//
//  Algorithm.swift
//  GasDynamics
//
//  Created by Арсений Токарев on 29.05.2022.
//

import Foundation

public class Algorithm1D: Algorithm {
    public typealias Mesh = Array<Cell>
    
    public let space   : Grid
    public var solution: [Double: Mesh] = [:]
    
    public override var data: Data? {
        var json: [String: Any] = [:]
        solution.sorted(by: { $0.key < $1.key }).forEach { key, value in
            guard Int(key)%35 == 0 else { return }
            let mesh = value.reduce(into: [String: String]()) {
                if let left = $1.left, let yL = left.y  {
                    $0[String(left.x)] = String(yL)
                }
                if let y = $1.middle.y {
                    $0[String($1.middle.x)] = String(y)
                }
                if let right = $1.right, let yR = right.y  {
                    $0[String(right.x)] = String(yR)
                }
            }
            json[String(key)] = mesh
        }
        return try? JSONSerialization.data(withJSONObject: json, options: .sortedKeys)
    }
    
    public final var interpolationData: Data? {
        var json: [String: Any] = [:]
        for t in stride(from: time.start, through: time.end, by: 35) {
            let mesh = stride(from: space.start, through: space.end, by: 0.1).reduce(into: [String: String]()) {
                guard let y = f(x: $1, t: t) else { return () }
                $0[String($1)] = String(y)
            }
            json[String(t)] = mesh
        }
        return try? JSONSerialization.data(withJSONObject: json, options: .sortedKeys)
    }
    
    public final var normC: Double {
        var max = Double.leastNormalMagnitude
        for t in time.nodes {
            guard let mesh = solution[t] else { continue }
            for cell in mesh {
                let middle = cell.middle.x
                let a = cell.left?.x ?? middle
                let b = cell.right?.x ?? middle
                let N = (b-a)/5
                for x in stride(from: a, through: b, by: N) {
                    guard let y = f(x: x, t: t) else { continue }
                    max = Swift.max(max, fabs(y))
                }
            }
        }
        return max
    }
//    public final var normL1: Double {
//        
//    }
    
    public required init(a: Double, b: Double, h: Double, tau: Double, deadline: Double) {
        self.space = Grid(start: a, end: b, step: h)
        super.init(tau: tau, deadline: deadline)
    }
    public required init(tau: Double = 1.0, deadline: Double) {
        self.space = .empty
        super.init(tau: tau, deadline: deadline)
    }
    
    public final func f(x: Double, t: Double) -> Double? {
        guard x >= space.start && x <= space.end else { return nil }
        guard let mesh = solution[t],
              let cell = mesh.first(where: {
                    let lBound = ($0.left?.x ?? $0.middle.x) <= x
                    let rBound = ($0.right?.x ?? $0.middle.x) >= x
                    return lBound && rBound
              }),
              let xi = xi(for: x, in: cell),
              let delta = delta(for: cell),
              let sixth = sixth(for: cell),
              let yL = cell.left?.y ?? cell.middle.y
        else { return nil }
        return yL+xi*(delta+sixth*(1-xi))
    }
    
    public final func adjust(mesh: inout Mesh) {
        for (i, cell) in mesh.enumerated() {
            guard let left = cell.left?.y,
                  let middle = cell.middle.y,
                  let right = cell.right?.y
            else { continue }
            guard (middle-left)*(right-middle) > 0 else {
                mesh[i].left?.y = middle
                mesh[i].right?.y = middle
                continue
            }
            guard let delta = delta(for: cell),
                  let six = sixth(for: cell)
            else { continue }
            let deltaSix = delta*six
            let deltaSq = pow(delta, 2)
            if deltaSix > deltaSq {
                mesh[i].left?.y = 3*middle-2*right
            }
            if deltaSix < -deltaSq {
                mesh[i].right?.y = 3*middle-2*left
            }
        }
    }
}

//MARK: - Grid
extension Algorithm {
    public final class Grid: Equatable {
        //MARK: Stepper
        /// - parameter:  order of node (Int)
        /// - result: value of node (Double)
        public typealias StepProvider = (Int) -> Double
        
        public let start : Double
        public let end   : Double
        public let steps : Int
        
        public final var step: Double {
            return (end-start)/Double(steps+1)
        }
        
        private var _nodes: [Double]?
        public final var nodes: [Double] {
            guard let cached = _nodes else {
                var values: [Double] = []
                values.reserveCapacity(steps)
                for x in stride(from: start, through: end, by: step) {
                    values.append(x)
                }
                _nodes = values
                return values
            }
            return cached
        }
        
        public init(start: Double, end: Double, steps: Int) {
            self.start = start
            self.end = end
            self.steps = steps
        }
        public init(start: Double, end: Double, step: Double) {
            self.start = start
            self.end = end
            self.steps = max(0, Int((end-start)/step) - 1)
        }
        
        public final func node(for i: Int) -> Double? {
            let node = Double(i)*step
            guard node >= start && node <= end else { return nil }
            return node
        }
        
        public static func == (lhs: Algorithm.Grid, rhs: Algorithm.Grid) -> Bool {
            return lhs.start == rhs.start && lhs.end == rhs.end && lhs.steps == rhs.steps
        }
    }
}

extension Algorithm1D.Grid {
    public static let empty = Algorithm1D.Grid(start: 0.0, end: 0.0, step: 0.0)
}

//MARK: - X -> Y Mapper
extension Algorithm1D {
    public struct Cell: Hashable {
        public let id = UUID()
        public var left  : Pair?
        public var middle: Pair
        public var right : Pair?
        
        public init(left: Pair? = nil, middle: Pair, right: Pair?) {
            self.left = left
            self.middle = middle
            self.right = right
        }
        
        public static func == (lhs: Algorithm1D.Cell, rhs: Algorithm1D.Cell) -> Bool {
            return lhs.id == rhs.id
        }
        
        public func hash(into hasher: inout Hasher) {
            hasher.combine(id)
        }
        
        public struct Pair {
            let x: Double
            var y: Double?
        }
    }
}

//MARK: - Helpers
extension Algorithm1D {
    public final func xi(for x: Double, in cell: Cell) -> Double? {
        let xL = cell.left?.x ?? cell.middle.x
        return (x-xL)/space.step
    }
    public final func delta(for cell: Cell) -> Double? {
        guard let left = cell.left?.y ?? cell.middle.y,
              let right = cell.right?.y ?? cell.middle.y
        else { return nil }
        return right-left
    }
    public final func sixth(for cell: Cell) -> Double? {
        guard let middle = cell.middle.y else { return nil }
        let left = cell.left?.y ?? middle
        let right = cell.right?.y ?? middle
        return 6*(middle-0.5*(left+right))
    }
}
