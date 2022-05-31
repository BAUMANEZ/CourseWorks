//
//  Advection.swift
//  GasDynamics
//
//  Created by Арсений Токарев on 31.05.2022.
//

import Foundation

public final class Advection1D: Algorithm1D {
    public let profile: Profile
    
    public let c    : Double
    public let l1   : Double
    public let l11  : Double
    public let l12  : Double
    public let l22  : Double
    public let l2   : Double
    public let L    : Double
    public let gamma: Double
    
    public init(
        profile : Profile,
        c       : Double = 1.0,
        l1      : Double = 10.0,
        l11     : Double = 50/3,
        l12     : Double = 20.0,
        l22     : Double = 70/3,
        l2      : Double = 30.0,
        L       : Double = 200.0,
        h       : Double = 1.0,
        gamma   : Double = 1.0,
        deadline: Double = 400.0
    ) {
        self.profile = profile
        self.c   = c
        self.l1  = l1
        self.l11 = l11
        self.l12 = l12
        self.l22 = l22
        self.l2  = l2
        self.L   = L
        self.gamma = gamma
        super.init(a: l1, b: L, h: h, tau: gamma*h/fabs(c), deadline: deadline)
        solve()
    }
    public required init(a: Double, b: Double, h: Double, tau: Double, deadline: Double) {
        self.profile = .rectangle
        self.c = 1.0
        self.l1  = a
        self.l11 = (b-a)/3
        self.l12 = (b-a)/3
        self.l22 = (b-a)/3
        self.l2  = (b-a)/2
        self.L   = b
        self.gamma = 1.0
        super.init(a: a, b: b, h: h, tau: tau, deadline: deadline)
    }
    public required init(tau: Double = 1.0, deadline: Double) {
        self.profile = .rectangle
        self.c = 0.0
        self.l1  = 0.0
        self.l11 = 0.0
        self.l12 = 0.0
        self.l22 = 0.0
        self.l2  = 0.0
        self.L   = 0.0
        self.gamma = 0.0
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
            var mesh: Mesh = [self.cell(for: x0, explicit: .zero)]
            for (index, cell) in prevMesh[1...].enumerated() {
                guard prevMesh.indices.contains(index),
                      let y = cell.middle.y,
                      let yM = prevMesh[index].middle.y
                else { continue }
                mesh.append(self.cell(for: cell.middle.x, explicit: (1-gamma)*y + gamma*yM))
            }
            solution[t] = mesh
            previousT = t
        }
        save(file: "advection", data: data)
    }
    private final func cell(for x: Double, explicit y: Double? = nil) -> Cell {
        let left = x-0.5 >= space.start ? Cell.Pair(x: x-0.5, y: nil) : nil
        let middle = Cell.Pair(x: x, y: y ?? profile.value(for: x, with: self))
        let right = x+0.5 <= space.end ? Cell.Pair(x: x+0.5, y: nil) : nil
        return Cell(left: left, middle: middle, right: right)
    }
}

//MARK: - Profiles
extension Advection1D {
    public enum Profile {
        case leftTriangle
        case rightTriangle
        case rectangle
        case tooth
        case M
        case cos
        
        public func value(for x: Double, with algorithm: Advection1D) -> Double {
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
