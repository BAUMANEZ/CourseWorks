//
//  Algorithm.swift
//  GasDynamics
//
//  Created by Арсений Токарев on 29.05.2022.
//

import Foundation

public class Algorithm {
    public let time: Grid
    
    public var data: Data? {
        return nil
    }
    
    public required init(tau: Double = 1.0, deadline: Double) {
        self.time = Grid(start: 0, end: deadline, step: tau)
        solve()
    }
    
    public func solve() {}
    
    public func save(file name: String) {
        let manager = FileManager.default
        let folder = manager.homeDirectoryForCurrentUser.appendingPathComponent("Desktop", isDirectory: true).appendingPathComponent("CourseWorks", isDirectory: true).appendingPathComponent("Gas Dynamics", isDirectory: true).appendingPathComponent("Code", isDirectory: true).appendingPathComponent("data", isDirectory: true)
        let file = folder.appendingPathComponent(name).appendingPathExtension("json")
        guard let _ = try? manager.createDirectory(at: folder, withIntermediateDirectories: true) else { return }
        manager.createFile(atPath: file.path, contents: data)
    }
}

public class Algorithm1D: Algorithm {
    public typealias Mesh = Array<Cell>
    
    public let space  : Grid
    public let profile: Profile
    
    public var solution: [Double: Mesh] = [:]
    
    public let a  : Double
    public let l1 : Double
    public let l11: Double
    public let l12: Double
    public let l22: Double
    public let l2 : Double
    public let L  : Double
    
    public var gamma: Double {
        return a*time.step/space.step
    }
    
    public override var data: Data? {
        var json: [String: Any] = [:]
        solution.sorted(by: { $0.key < $1.key }).forEach { key, value in
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
    
    public init(
        profile : Profile,
        a       : Double = 1.0,
        l1      : Double = 10.0,
        l11     : Double = 50/3,
        l12     : Double = 20.0,
        l22     : Double = 70/3,
        l2      : Double = 30,
        L       : Double = 150,
        h       : Double = 1,
        deadline: Double = 400
    ) {
        self.profile = profile
        self.a = a
        self.l1  = l1
        self.l11 = l11
        self.l12 = l12
        self.l22 = l22
        self.l2  = l2
        self.L   = L
        let nodes = Grid(start: l1, end: L, step: h)
        self.space = nodes
        super.init(tau: h / fabs(a), deadline: deadline)
    }
    
    public required init(tau: Double = 1.0, deadline: Double) {
        self.profile = .rectangle
        self.a = 0.0
        self.l1  = 0.0
        self.l11 = 0.0
        self.l12 = 0.0
        self.l22 = 0.0
        self.l2  = 0.0
        self.L   = 0.0
        self.space = Grid(start: 0.0, end: 0.0, step: 0.0)
        super.init(tau: tau, deadline: deadline)
    }

    public override func solve() {
        guard time.steps > 1,
              let t0 = time.nodes.first
        else { return }
        solution[t0] = space.nodes.reduce(into: Mesh()) { $0.append(self.cell(for: $1)) }
        let gamma = self.gamma
        var previousT = t0
        for t in time.nodes[1...] {
            guard
                let prevMesh = solution[previousT],
                let x0 = prevMesh.first?.middle.x
            else { continue }
            var mesh: Mesh = [self.cell(for: x0, exlpicit: .zero)]
            for (index, cell) in prevMesh[1...].enumerated() {
                guard prevMesh.indices.contains(index-1),
                      let y = cell.middle.y,
                      let yM = prevMesh[index-1].middle.y
                else { continue }
                mesh.append(self.cell(for: cell.middle.x, exlpicit: (1-gamma)*y + gamma*yM))
            }
            solution[t] = mesh
            previousT = t
        }
        save(file: "advection")
    }
    
    private func cell(for x: Double, exlpicit y: Double? = nil) -> Cell {
        let left = x-0.5 >= l1 ? Cell.Pair(x: x-0.5, y: nil) : nil
        let middle = Cell.Pair(x: x, y: y ?? profile.value(for: x, with: self))
        let right = x+0.5 <= L ? Cell.Pair(x: x+0.5, y: nil) : nil
        return Cell(left: left, middle: middle, right: right)
    }
}

extension Algorithm {
    public class Grid: Equatable {
        //MARK: Stepper
        /// - parameter:  order of node (Int)
        /// - result: value of node (Double)
        public typealias StepProvider = (Int) -> Double
        
        public let start : Double
        public let end   : Double
        public let steps : Int
        
        public var step: Double {
            return (end-start)/Double(steps+1)
        }
        
        private var _nodes: [Double]?
        public var nodes: [Double] {
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
        
        public func node(for i: Int) -> Double? {
            let node = Double(i)*step
            guard node >= start && node <= end else { return nil }
            return node
        }
        
        public static func == (lhs: Algorithm.Grid, rhs: Algorithm.Grid) -> Bool {
            return lhs.start == rhs.start && lhs.end == rhs.end && lhs.steps == rhs.steps
        }
    }
}

extension Algorithm1D {
    public struct Cell: Hashable {
        public var left  : Pair?
        public var middle: Pair
        public var right : Pair?
        
        public init(left: Pair? = nil, middle: Pair, right: Pair?) {
            self.left = left
            self.middle = middle
            self.right = right
        }
        
        public static func == (lhs: Algorithm1D.Cell, rhs: Algorithm1D.Cell) -> Bool {
            let left = lhs.left?.x == rhs.left?.x && lhs.left?.y == rhs.left?.y
            let middle = lhs.middle.x == rhs.middle.x && lhs.middle.y == rhs.middle.y
            let right = lhs.right?.x == rhs.right?.x && lhs.right?.y == rhs.right?.y
            return left && middle && right
        }
        
        public func hash(into hasher: inout Hasher) {
            hasher.combine(left?.x)
            hasher.combine(left?.y)
            hasher.combine(middle.x)
            hasher.combine(middle.y)
            hasher.combine(right?.x)
            hasher.combine(right?.y)
        }
        
        public struct Pair {
            var x: Double
            var y: Double?
        }
    }
}

extension Algorithm1D {
    public func delta(for cell: Cell) -> Double? {
        guard let left = cell.left?.y,
              let right = cell.right?.y
        else { return nil }
        return right-left
    }
    
    public func sixth(for cell: Cell) -> Double? {
        guard let left = cell.left?.y,
              let middle = cell.middle.y,
              let right = cell.right?.y
        else { return nil }
        return 6*(middle-0.5*(left+right))
    }
}


extension Algorithm1D {
    public enum Profile {
        case leftTriangle
        case rightTriangle
        case rectangle
        case tooth
        case M
        case cos
        
        public func value(for x: Double, with algorithm: Algorithm1D) -> Double {
            guard x >= algorithm.l1 && x <= algorithm.l2 else { return .zero }
            switch self {
            case .leftTriangle:
                return (x-algorithm.l1)/(algorithm.l2-algorithm.l1)
            case .rightTriangle:
                return (algorithm.l2-x)/(algorithm.l2-algorithm.l1)
            case .rectangle:
                return 1.0
            case .cos:
                return 0.5-0.5*Foundation.cos(2*Double.pi/(algorithm.l2-algorithm.l1) * (x-algorithm.l1) )
            case .tooth:
                switch x {
                case algorithm.l1 ..< algorithm.l11:
                    return -2*(x-algorithm.l1)/(3*(algorithm.l11-algorithm.l1)) + 1
                case algorithm.l11 ... algorithm.l22:
                    return 1/3
                default:
                    return 2*(x-algorithm.l2)/(3*(algorithm.l2-algorithm.l22)) + 1
                }
            case .M:
                switch x {
                case algorithm.l1 ..< algorithm.l12:
                    return -2*(x-algorithm.l1)/(3*(algorithm.l12-algorithm.l1)) + 1
                default:
                    return 2*(x-algorithm.l2)/(3*(algorithm.l2-algorithm.l12)) + 1
                }
            }
        }
    }
}
