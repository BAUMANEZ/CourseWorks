//
//  Algorithm.swift
//  GasDynamics
//
//  Created by Арсений Токарев on 29.05.2022.
//

import Foundation

public class Algorithm {
    public let time: Grid
    
    public required init(tau: Double = 1.0, deadline: Double) {
        self.time = Grid(start: 0, end: deadline, step: tau)
    }
}

public class Algorithm1D: Algorithm {
    public typealias Mesh = [Boundary: Boundary]
    
    public let space  : Grid
    public let profile: Profile
    public let segments: [Boundary]
    
    // Time: Mesh
    public private(set) var solution: [Double: Mesh] = [:]
    
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
        self.segments = nodes.nodes.map {
            let left = $0-0.5 >= l1 ? $0-0.5 : nil
            let right = $0+0.5 <= L ? $0+0.5 : nil
            return Boundary(left: left, middle: $0, right: right)
        }
        super.init(tau: h / fabs(a), deadline: deadline)
        solve()
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
        self.segments = []
        super.init(tau: tau, deadline: deadline)
    }
    
    private func solve() {
        guard time.steps > 1,
              let t0 = time.nodes.first,
              segments.count > 1,
              let x0 = segments.first
        else { return }
        solution[t0] = segments.reduce(into: Mesh()) {
            let value = profile.value(for: $1.middle, with: self)
            $0[$1] = Boundary(middle: value)
        }
        let gamma = self.gamma
        var previousT = t0
        for t in time.nodes[1...] {
            guard let prevMesh = solution[previousT] else { continue }
            var mesh: Mesh = [x0: Boundary(middle: .zero)]
            for (index, segment) in segments[1...].enumerated() {
                guard segments.indices.contains(index-1),
                      let y = prevMesh[segment],
                      let yMinus = prevMesh[segments[index-1]]
                else { continue }
                mesh[segment] = Boundary(middle: (1-gamma)*y.middle + gamma*yMinus.middle)
            }
            solution[t] = mesh
            previousT = t
        }
    }
}

extension Algorithm {
    public struct Grid: Equatable {
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
        
        public var nodes: [Double] {
            var values: [Double] = []
            values.reserveCapacity(steps)
            for x in stride(from: start, through: end, by: step) {
                values.append(x)
            }
            return values
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
    }
}

extension Algorithm {
    public struct Boundary: Hashable {
        public var left  : Double
        public var middle: Double
        public var right : Double
        
        public init(left: Double? = nil, middle: Double, right: Double? = nil) {
            self.left = left ?? middle
            self.middle = middle
            self.right = right ?? middle
        }
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
