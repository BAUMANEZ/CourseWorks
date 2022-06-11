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
        return 4.0
    }
    
    private final var plotStep: Int {
        return time.steps/5
    }
    
    private lazy var F: (Double) -> Double = { return self.c*$0 }
    
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
        detailed[time.start] = ([space.start-space.step]+space.nodes).reduce(into: DetailedMesh()) {
            let xL = Node(value: $1-space.halfed, side: .right)
            let xR = Node(value: $1+space.halfed, side: .left)
            let x  = Node(value: $1, side: .middle)
            let yL = self.u0(xL.value)
            let yR = self.u0(xR.value)
            let y  = (yL+yR)/2.0
            $0[xL] = yL
            $0[xR] = yR
            $0[x]  = y
        }
        guard time.steps > 1 else { return }
        time.nodes[1...].forEach{ solve(for: $0) }
        var solutions: [Time: Mesh] = [:]
        let tttt = 4
        for i in stride(from: 0, through: tttt*15, by: 15) {
            let t = time.nodes[i]
            guard detailed[t] != nil else { return () }
            solutions[t] = space.nodes.reduce(into: Mesh()) { mesh, node in
                let xL = Node(value: node-space.halfed, side: .right)
                let x  = Node(value: node, side: .middle)
                let xR = Node(value: node+space.halfed, side: .left)
                guard let y  = detailed[t]?[x],
                      let yL = detailed[t]?[xL],
                      let yR = detailed[t]?[xR]
                else { return }
                mesh[x.value] = y
                mesh[xL.value+10e-10] = yL
                mesh[xR.value-10e-10] = yR
//                let delta = delta(yL: yL, yR: yR)
//                let sixth = sixth(yL: yL, y: y, yR: yR)
//                for j in stride(from: xL.value, through: xR.value, by: (xR.value-xL.value)/substep) {
//                    let xi = (j-xL.value)/space.step
//                    let ttt = j == xL.value ? j+10e-11 : (j == xR.value ? j-10e-11 : j)
//                    mesh[ttt] = yL+xi*(delta+sixth*(1.0-xi))
//                }
            }
        }
        if let identifier = identifier {
            save(file: identifier, data: data(for: solutions), meta: metadata)
        }
    }
    
    public func solve(for t: Time) {
        detailed[t] = [:]
        let tP = t-time.step
        if t==24 {
            let a = 4
        }
        space.nodes.forEach { node in
            let x  = Node(value: node, side: .middle)
            let xL = Node(value: node-space.halfed, side: .right)
            let xR = Node(value: node+space.halfed, side: .left)
            guard let leftFlow = flow(t: tP, xB: xL),
                  let rightFlow = flow(t: tP, xB: xR),
                  let y = detailed[tP]?[x]
            else { return }
            let f = y-1.0/(space.step)*(rightFlow-leftFlow)
            detailed[t]?[x] = f
        }
    }
    
    private func flow(t: Time, xB: Node) -> Double? {
        let xL: Node
        let x : Node
        let xR: Node
        if c > 0 {
            xL = Node(value: xB.value-space.step, side: .right)
            x  = Node(value: xB.value-space.halfed, side: .middle)
            xR = Node(value: xB.value, side: .left)
        } else {
            xL = Node(value: xB.value, side: .right)
            x  = Node(value: xB.value+space.halfed, side: .middle)
            xR = Node(value: xB.value+space.step, side: .left)
        }
        guard let y = detailed[t]?[x],
              let yL = detailed[t]?[xL],
              let yR = detailed[t]?[xR]
        else { return nil }
        let h = c*time.step/space.step
        let delta = delta(yL: yL, yR: yR)
        let sixth = sixth(yL: yL, y: y, yR: yR)
        let avg = c > 0 ? yR-0.5*h*(delta-(1-2.0/3.0*h)*sixth) : yL+0.5*h*(delta-(1-2.0/3.0*h)*sixth)
        return F(avg)
    }
}

extension InterpolatedAdvection1D {
    public final func delta(t: Time, xL: Node, xR: Node) -> Double? {
        guard let yL = detailed[t]?[xL],
              let yR = detailed[t]?[xR]
        else { return nil }
        return delta(yL: yL, yR: yR)
    }
    public final func delta(yL: Double, yR: Double) -> Double {
        return yR-yL
    }
    public final func sixth(t: Time, xL: Node, x: Node, xR: Node) -> Double? {
        guard let y = detailed[t]?[x],
              let yL = detailed[t]?[xL],
              let yR = detailed[t]?[xR]
        else { return nil }
        return sixth(yL: yL, y: y, yR: yR)
    }
    public final func sixth(yL: Double, y: Double, yR: Double) -> Double {
        return 6.0*(y-0.5*(yR+yL))
    }
    public final func xi(x: Double, xL: Node) -> Double {
        return (x-xL.value)/space.step
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
    public final func Nf(x: Double, xL: Node, yL: Double, y: Double, yR: Double) -> Double {
        let xi = xi(x: x, xL: xL)
        let delta = delta(yL: yL, yR: yR)
        let sixth = sixth(yL: yL, y: y, yR: yR)
        return yL+xi*(delta+sixth*(1-xi))
    }
}
