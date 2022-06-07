//
//  Interpolated1D.swift
//  GasDynamics
//
//  Created by Арсений Токарев on 07.06.2022.
//

import Foundation

public class InterpolatedAdvection1D: Advection1D {
    public typealias DetailedMesh = [Node: Double]
    public var detailed: [Time: DetailedMesh] = [:]
    
    public var identifier: String? {
        return nil
    }
    
    public var substep: Double {
        return 20.0
    }
    
    public final var normC: Double {
        var result = -Double.greatestFiniteMagnitude
        time.nodes[1...].forEach { t in
            guard detailed[t] != nil else { return }
            space.nodes.forEach { node in
                let xL = Node(value: node-space.halfed, side: .right)
                let xM = Node(value: node, side: .middle)
                let xR = Node(value: node+space.halfed, side: .left)
                guard let y = detailed[t]?[xM],
                      let yL = detailed[t]?[xL],
                      let yR = detailed[t]?[xR]
                else { return }
                let delta = yR-yL
                let sixth = 6.0*(y-0.5*(yR+yL))
                for j in stride(from: xL.value, through: xR.value, by: (xR.value-xL.value)/200.0) {
                    let xi = (j-xL.value)/space.step
                    result = Swift.max(result, (u(j, t)-(yL+xi*(delta+sixth*(1.0-xi)))).magnitude)
                }
            }
        }
        return result
    }
    
    private final func Nf(in x: Double, xL: Double, yL: Double, y: Double, xR: Double, yR: Double) -> Double {
        let xi = (x-xL)/space.step
        let delta = yR-yL
        let sixth = 6.0*(y-0.5*(yR+yL))
        return yL+xi*(delta+sixth*(1.0-xi))
    }
    
    public final override func f(x: Double, t: Double) -> Double? {
        guard let _x = space.nodes.first(where: { $0-space.halfed <= x && $0+space.halfed >= x }) else { return nil }
        let xL = Node(value: _x-space.halfed, side: .right)
        let xM = Node(value: _x, side: .middle)
        let xR = Node(value: _x+space.halfed, side: .left)
        guard let y = detailed[t]?[xM],
              let yL = detailed[t]?[xL],
              let yR = detailed[t]?[xR]
        else { return nil }
        return Nf(in: x, xL: xL.value, yL: yL, y: y, xR: xR.value, yR: yR)
    }
    
    public final override func solve() {
        detailed[time.start] = space.nodes.reduce(into: DetailedMesh()) {
            $0[Node(value: $1, side: .middle)] = self.u0($1)
            let xL = $1-space.halfed
            if xL>=space.start { $0[Node(value: xL, side: .right)] = self.u0(xL) }
            let xR = $1+space.halfed
            if xR<=space.end { $0[Node(value: xR, side: .left)] = self.u0(xR) }
        }
        guard time.steps > 1 else { return }
        time.nodes[1...].forEach{ solve(for: $0) }
        adjust()
        
        let solutions = time.nodes[1...].reduce(into: [Time: Mesh]()) { solutions, t in
            guard detailed[t] != nil else { return () }
            solutions[t] = space.nodes.reduce(into: Mesh()) { mesh, node in
                let xL = Node(value: node-space.halfed, side: .right)
                let xM = Node(value: node, side: .middle)
                let xR = Node(value: node+space.halfed, side: .left)
                guard let y = detailed[t]?[xM],
                      let yL = detailed[t]?[xL],
                      let yR = detailed[t]?[xR]
                else { return }
                let delta = yR-yL
                let sixth = 6.0*(y-0.5*(yR+yL))
                for j in stride(from: xL.value, through: xR.value, by: (xR.value-xL.value)/substep) {
                    let xi = (j-xL.value)/space.step
                    mesh[j] = yL+xi*(delta+sixth*(1.0-xi))
                }
            }
        }
        if let identifier = identifier {
            save(file: identifier, data: data(for: solutions))
        }
//        let solutions = time.nodes[1...].reduce(into: [Time: Mesh]()) { solutions, t in
//            var mesh: Mesh = [:]
//            for x in stride(from: space.start+2.0*space.step, through: space.end, by: space.step/substep) {
//                mesh[x] = f(x: x, t: t)
//            }
//            solutions[t] = mesh
//        }
    }
    
    public func solve(for t: Time) {
        detailed[t] = space.nodes.reduce(into: DetailedMesh()) {
            $0[Node(value: $1, side: .middle)] = drift(for: t, in: $1)
        }
    }
    
    private final func adjust() {
        for t in time.nodes {
            guard detailed[t] != nil else { continue }
            for node in space.nodes {
                let xL = Node(value: node-space.halfed, side: .right)
                let x = Node(value: node, side: .middle)
                let xR = Node(value: node+space.halfed, side: .left)
                guard let yL = detailed[t]?[xL],
                      let y = detailed[t]?[x],
                      let yR = detailed[t]?[xR]
                else { continue }
                guard (yR-y)*(y-yL) > 0 else {
                    detailed[t]?[xL] = y
                    detailed[t]?[xR] = y
                    continue
                }
                let delta = yR-yL
                let sixth = 6.0*(y-0.5*(yL+yR))
                let deltaSix = delta*sixth
                let deltaSq = delta*delta
                if deltaSix > deltaSq {
                    detailed[t]?[xL] = 3.0*y-2.0*yR
                }
                if deltaSix < -deltaSq {
                    detailed[t]?[xR] = 3.0*y-2.0*yL
                }
            }
        }
    }
}
