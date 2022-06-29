//
//  Interpolated1D.swift
//  GasDynamics
//
//  Created by Арсений Токарев on 07.06.2022.
//

import Foundation

public class InterpolatedAdvection1D: Advection1D {
    public typealias DetailedMesh = [Node: Double]
    public var detailed: [Int: DetailedMesh] = [:]
    
    public var identifier: String? {
        return nil
    }
    
    private final var substep: Double {
        return 4.0
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
                "C"  : String(normC),
                "L"  : nil,//String(normL),
                "L2" : String(normL2),
                "W2" : nil
            ]
        ]
        return try? JSONSerialization.data(withJSONObject: json, options: [.prettyPrinted, .sortedKeys])
    }
    
    public final var normC: Double {
        var result = -Double.greatestFiniteMagnitude
        time.range(starting: 1).forEach { j in
            guard detailed[j] != nil else { return }
            space.nodes().forEach { node in
                let xL = Node(value: node-space.halfed, side: .right)
                let xM = Node(value: node, side: .middle)
                let xR = Node(value: node+space.halfed, side: .left)
                guard let y = detailed[j]?[xM],
                      let yL = detailed[j]?[xL],
                      let yR = detailed[j]?[xR]
                else { return }
                let delta = yR-yL
                let sixth = 6.0*(y-0.5*(yR+yL))
                for k in stride(from: xL.value, to: xR.value, by: (xR.value-xL.value)/substep) {
                    let xi = (k-xL.value)/space.step
                    result = Swift.max(result, (u(k, time.node(for: j))-(yL+xi*(delta+sixth*(1.0-xi)))).magnitude)
                }
            }
        }
        return result
    }
    
    public final var normL: Double {
        time.range(starting: 1).reduce(into: Double()) { global, j in
            guard detailed[j] != nil else { return () }
            global += space.nodes().reduce(into: Double()) { local, node in
                let jP = j-1
                let xL = Node(value: node-space.halfed, side: .right)
                let x  = Node(value: node, side: .middle)
                let xR = Node(value: node+space.halfed, side: .left)
                guard let y = detailed[j]?[x],
                      let yL = detailed[j]?[xL],
                      let yR = detailed[j]?[xR],
                      let yP = detailed[jP]?[x],
                      let yPL = detailed[jP]?[xL],
                      let yPR = detailed[jP]?[xR]
                else { return }
                let delta = yR-yL
                let deltaP = yPR-yPL
                let sixth = 6.0*(y-0.5*(yR+yL))
                let sixthP = 6.0*(yP-0.5*(yPR+yPL))
                for k in stride(from: xL.value, through: xR.value, by: (xR.value-xL.value)/substep) {
                    let xi = (k-xL.value)/space.step
                    let z = Swift.abs(u(k, time.node(for: j))-(yL+xi*(delta+sixth*(1.0-xi))))
                    let zP = Swift.abs(u(k, time.node(for: jP))-(yPL+xi*(deltaP+sixthP*(1.0-xi))))
                    local += space.step*time.step*(zP-z)
                }
            }
        }
    }
    public final var normL2: Double {
        return sqrt(time.range(starting: 1).reduce(into: Double()) { global, j in
            guard detailed[j] != nil else { return () }
            global += space.nodes().reduce(into: Double()) { local, node in
                let jP = j-1
                let xL = Node(value: node-space.halfed, side: .right)
                let x  = Node(value: node, side: .middle)
                let xR = Node(value: node+space.halfed, side: .left)
                guard let y = detailed[j]?[x],
                      let yL = detailed[j]?[xL],
                      let yR = detailed[j]?[xR],
                      let yP = detailed[jP]?[x],
                      let yPL = detailed[jP]?[xL],
                      let yPR = detailed[jP]?[xR]
                else { return }
                let delta = yR-yL
                let deltaP = yPR-yPL
                let sixth = 6.0*(y-0.5*(yR+yL))
                let sixthP = 6.0*(yP-0.5*(yPR+yPL))
                for k in stride(from: xL.value, through: xR.value, by: (xR.value-xL.value)/substep) {
                    let xi = (k-xL.value)/space.step
                    let z = pow((u(k, time.node(for: j))-(yL+xi*(delta+sixth*(1.0-xi)))), 2)
                    let zP = pow(Swift.abs(u(k, time.node(for: jP))-(yPL+xi*(deltaP+sixthP*(1.0-xi)))), 2)
                    local += space.step*time.step*(zP-z)
                }
            }
        })
    }
    public final var normW2: Double {
        return .zero
    }
    
    public final override func solve() {
        detailed[0] = space.nodes().reduce(into: DetailedMesh()) {
            let xL = Node(value: $1-space.halfed, side: .right)
            let x  = Node(value: $1, side: .middle)
            let xR = Node(value: $1+space.halfed, side: .left)
            let yL = u0(xL.value)
            let yR = u0(xR.value)
            let y  = (yL+yR)/2.0
            $0[xL] = yL
            $0[x]  = y
            $0[xR] = yR
        }
        guard time.steps > 1 else { return }
        time.range(starting: 1).forEach{ solve(for: $0) }
        var solutions: [Time: Mesh] = [:]
        for j in Swift.stride(from: 0, through: time.steps, by: plotStep) {
//            let limit = 40
//            let stride = 1
//        for j in Swift.stride(from: 0, through: limit*stride, by: stride) {
            guard detailed[j] != nil else { continue }
            let t = time.node(for: j)
            solutions[t] = space.nodes().reduce(into: Mesh()) { mesh, node in
                let xL = Node(value: node-space.halfed, side: .right)
                let x  = Node(value: node, side: .middle)
                let xR = Node(value: node+space.halfed, side: .left)
                guard let yL = detailed[j]?[xL],
                      let y  = detailed[j]?[x],
                      let yR = detailed[j]?[xR]
                else { return }
                mesh[x.value] = y
                let delta = delta(yL: yL, yR: yR)
                let sixth = sixth(yL: yL, y: y, yR: yR)
                for k in Swift.stride(from: xL.value, through: xR.value, by: (xR.value-xL.value)/substep) {
                    let xi = (k-xL.value)/space.step
                    let ttt = k == xL.value ? k+10e-11 : (k == xR.value ? k-10e-11 : k)
                    mesh[ttt] = yL+xi*(delta+sixth*(1.0-xi))
                }
            }
        }
        if let identifier = identifier {
            save(file: identifier, data: data(for: solutions), meta: metadata)
        }
    }
    
    public func solve(for j: Int) {
        detailed[j] = [:]
        let jP = j-1
        space.nodes().forEach { node in
            let x  = Node(value: node, side: .middle)
            let xL = Node(value: node-space.halfed, side: .right)
            let xR = Node(value: node+space.halfed, side: .left)
            guard let leftFlow = flow(j: jP, xB: xL),
                  let rightFlow = flow(j: jP, xB: xR),
                  let y = detailed[jP]?[x]
            else { return }
            let f = y-gamma*(rightFlow-leftFlow)
            detailed[j]?[x] = f
        }
    }
    
    private func flow(j: Int, xB: Node) -> Double? {
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
        guard let y  = detailed[j]?[x],
              let yL = detailed[j]?[xL],
              let yR = detailed[j]?[xR]
        else { return .zero }
//        let delta = delta(yL: yL, yR: yR)
//        let sixth = sixth(yL: yL, y: y, yR: yR)
//        let avg = c > 0 ? yR-0.5*gamma*(delta-sixth*(1.0-2.0/3.0*gamma)) : yL+0.5*gamma*(delta-sixth*(1-2.0/3.0*gamma))
//        return c*avg
//
        let boundary = c > 0 ? xR.value : xL.value
        let f1 = c*(Nf(x: boundary-c*time.step, xL: xL, yL: yL, y: y, yR: yR))
        let f2 = c*(Nf(x: boundary-0.5*c*time.step, xL: xL, yL: yL, y: y, yR: yR))
        let f3 = c*(Nf(x: boundary, xL: xL, yL: yL, y: y, yR: yR))
        let flow = 1.0/6.0*(f1+4.0*f2+f3)
        return flow
    }
}

extension InterpolatedAdvection1D {
    public final func delta(j: Int, xL: Node, xR: Node) -> Double? {
        guard let yL = detailed[j]?[xL],
              let yR = detailed[j]?[xR]
        else { return nil }
        return delta(yL: yL, yR: yR)
    }
    public final func delta(yL: Double, yR: Double) -> Double {
        return yR-yL
    }
    public final func sixth(j: Int, xL: Node, x: Node, xR: Node) -> Double? {
        guard let y = detailed[j]?[x],
              let yL = detailed[j]?[xL],
              let yR = detailed[j]?[xR]
        else { return nil }
        return sixth(yL: yL, y: y, yR: yR)
    }
    public final func sixth(yL: Double, y: Double, yR: Double) -> Double {
        return 6.0*(y-0.5*(yR+yL))
    }
    public final func xi(x: Double, xL: Node) -> Double {
        return (x-xL.value)/space.step
    }
    public final func Nf(j: Int, x: Double, xL: Node, xM: Node, xR: Node) -> Double? {
        guard let y = detailed[j]?[xM],
              let yL = detailed[j]?[xL],
              let yR = detailed[j]?[xR]
        else { return .zero }
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
