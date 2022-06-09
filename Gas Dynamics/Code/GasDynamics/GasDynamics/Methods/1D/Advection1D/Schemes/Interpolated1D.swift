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
    
    private final var substep: Double {
        return 8.0
    }
    
    private final var plotStep: Int {
        return time.steps/5
    }
    
    private final var metadata: Data? {
        let json: [String: Any] = [
            "type"   : "Advection",
            "scheme" : identifier ?? "none",
            "profile": profile ?? "unknown",
            "grid"   : [
                "time" : time.json,
                "space": space.json
            ],
            "speed": c,
            "norms"  : [
                "C"  : nil,//String(normC),
                "L"  : nil,//String(normL),
                "L2" : String(normL2),
                "W2" : nil
            ]
        ]
        return try? JSONSerialization.data(withJSONObject: json, options: [.prettyPrinted, .sortedKeys])
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
        return .zero
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
                    let z = pow((u(j, t)-(yL+xi*(delta+sixth*(1.0-xi)))), 2)
                    let zP = pow(Swift.abs(u(j, tP)-(yPL+xi*(deltaP+sixthP*(1.0-xi)))), 2)
                    local += space.step*time.step*(zP-z)
                }
            }
        })
    }
    public final var normW2: Double {
        return .zero
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
        var solutions: [Time: Mesh] = [:]
        for t in stride(from: time.start+time.step, through: time.end, by: Double(plotStep)) {
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
            save(file: identifier, data: data(for: solutions), meta: metadata)
        }
    }
    
    public func solve(for t: Time) {
        let drift = c*t
        detailed[t] = [:]
        space.nodes[1...].forEach { node in
            let x  = Node(value: node, side: .middle)
            let xL = node-space.halfed
            let xR = node+space.halfed
            let fa = u0(xL-drift)
            let fb = u0(xR-drift)
            let fh = u0((xL+xR)/2.0-drift)
            let y = 1.0/6.0*(fa+4.0*fh+fb)
            detailed[t]?[x] = y
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
