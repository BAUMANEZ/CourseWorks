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
        return 100.0
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
                for j in stride(from: xL.value, to: xR.value, by: (xR.value-xL.value)/substep) {
                    let xi = (j-xL.value)/space.step
                    result = Swift.max(result, (u(j, t)-(yL+xi*(delta+sixth*(1.0-xi)))).magnitude)
                }
            }
        }
        return result
    }
    
    public final var normL: Double {
        return time.nodes[0..<time.nodes.count-1].reduce(into: Double()) { global, t in
            guard detailed[t] != nil else { return () }
            global += space.nodes.reduce(into: Double()) { local, node in
                let tP = t+time.step
                let xL = Node(value: node-space.halfed, side: .right)
                let x  = Node(value: node, side: .middle)
                let xR = Node(value: node+space.halfed, side: .left)
                guard let y = detailed[t]?[x],
                      let yL = detailed[t]?[xL],
                      let yR = detailed[t]?[xR],
                      let yP = detailed[tP]?[x],
                      let yPL = detailed[tP]?[xL],
                      let yPR = detailed[tP]?[xR]
                else { return }
                let delta = yR-yL
                let deltaP = yPR-yPL
                let sixth = 6.0*(y-0.5*(yR+yL))
                let sixthP = 6.0*(yP-0.5*(yPR+yPL))
                for j in stride(from: xL.value, through: xR.value, by: (xR.value-xL.value)/substep) {
                    let xi = (j-xL.value)/space.step
                    let f = Swift.abs(u(j, t)-(yL+xi*(delta+sixth*(1.0-xi))))
                    let fP = Swift.abs(u(j, tP)-(yPL+xi*(deltaP+sixthP*(1.0-xi))))
                    local += space.step*(fP-f)
                }
            }
        }
    }
    public final var normL2: Double {
        return sqrt(time.nodes[0..<time.nodes.count-1].reduce(into: Double()) { global, t in
            guard detailed[t] != nil else { return () }
            global += space.nodes.reduce(into: Double()) { local, node in
                let tP = t+time.step
                let xL = Node(value: node-space.halfed, side: .right)
                let x  = Node(value: node, side: .middle)
                let xR = Node(value: node+space.halfed, side: .left)
                guard let y = detailed[t]?[x],
                      let yL = detailed[t]?[xL],
                      let yR = detailed[t]?[xR],
                      let yP = detailed[tP]?[x],
                      let yPL = detailed[tP]?[xL],
                      let yPR = detailed[tP]?[xR]
                else { return }
                let delta = yR-yL
                let deltaP = yPR-yPL
                let sixth = 6.0*(y-0.5*(yR+yL))
                let sixthP = 6.0*(yP-0.5*(yPR+yPL))
                for j in stride(from: xL.value, through: xR.value, by: (xR.value-xL.value)/substep) {
                    let xi = (j-xL.value)/space.step
                    let f = Swift.abs(u(j, t)-(yL+xi*(delta+sixth*(1.0-xi))))
                    let fP = Swift.abs(u(j, tP)-(yPL+xi*(deltaP+sixthP*(1.0-xi))))
                    local += pow(space.step*(fP-f), 2)
                }
            }
        })
    }
    public final var normW2: Double {
        return sqrt(time.nodes[0 ..< time.nodes.count-1].reduce(into: Double()) { global, t in
            global += space.nodes.reduce(into: Double()) { local, node in
                let x = Node(value: node, side: .middle)
                guard let y = detailed[t]?[x],
                      let yP = detailed[t+time.step]?[x]
                else { return () }
                local += pow(space.step*(abs(yP-u(node, t)-abs(y-u(node, t+time.step)))), 2)
            }
        })
    }
    
    public final override func solve() {
        detailed[time.start] = space.nodes.reduce(into: DetailedMesh()) {
            $0[Node(value: $1, side: .middle)] = self.u0($1)
            let xL = $1-space.halfed
            $0[Node(value: xL, side: .right)] = self.u0(xL)
            let xR = $1+space.halfed
            $0[Node(value: xR, side: .left)] = self.u0(xR)
        }
        guard time.steps > 1 else { return }
        time.nodes[1...].forEach{ solve(for: $0) }
        let solutions = time.nodes[1...].reduce(into: [Time: Mesh]()) { solutions, t in
            guard detailed[t] != nil else { return () }
            solutions[t] = space.nodes.reduce(into: Mesh()) { mesh, node in
                let xL = Node(value: node-space.halfed, side: .right)
                let x  = Node(value: node, side: .middle)
                let xR = Node(value: node+space.halfed, side: .left)
                guard let y  = detailed[t]?[x],
                      let yL = detailed[t]?[xL],
                      let yR = detailed[t]?[xR]
                else { return }
                let delta = yR-yL
                let sixth = 6.0*(y-0.5*(yR+yL))
                for j in stride(from: xL.value, through: xR.value, by: (xR.value-xL.value)/substep) {
                    let xi = (j-xL.value)/space.step
                    let ttt = j == xL.value ? j+10e-11 : (j == xR.value ? j-10e-11 : j)
                    mesh[ttt] = yL+xi*(delta+sixth*(1.0-xi))
                }
            }
        }
        if let identifier = identifier {
            save(file: identifier, data: data(for: solutions))
        }
    }
    
    public func solve(for t: Time) {
        detailed[t] = space.nodes.reduce(into: DetailedMesh()) {
            $0[Node(value: $1, side: .middle)] = drift(for: t, in: $1)
        }
    }
}

extension InterpolatedAdvection1D {
    public final func delta(t: Time, xL: Node, xR: Node) -> Double? {
        guard let yL = detailed[t]?[xL],
              let yR = detailed[t]?[xR]
        else { return nil }
        return yR-yL
    }
    public final func sixth(t: Time, xL: Node, x: Node, xR: Node) -> Double? {
        guard let y = detailed[t]?[x],
              let yL = detailed[t]?[xL],
              let yR = detailed[t]?[xR]
        else { return nil }
        return 6.0*(y-0.5*(yR+yL))
    }
    public final func Nf(t: Time, x: Double, xL: Node, xM: Node, xR: Node) -> Double? {
        guard let y = detailed[t]?[xM],
              let yL = detailed[t]?[xL],
              let yR = detailed[t]?[xR]
        else { return nil }
        let xi = (x-xL.value)/space.step
        let delta = yR-yL
        let sixth = 6.0*(y-0.5*(yR+yL))
        return yL+xi*(delta+sixth*(1.0-xi))
    }
}
