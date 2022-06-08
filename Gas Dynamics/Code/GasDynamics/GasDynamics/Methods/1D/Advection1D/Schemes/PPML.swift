//
//  PPML.swift
//  GasDynamics
//
//  Created by Арсений Токарев on 02.06.2022.
//

import Foundation

extension Advection1D {
    public final class PPML: InterpolatedAdvection1D {
        public override var identifier: String? {
            return "advectionPPML"
        }
        
        public override func solve(for t: Time) {
            super.solve(for: t)
            let tP = t-time.step
            guard detailed[t] != nil, detailed[tP] != nil else { return }
            space.nodes.forEach { node in
                let x  = Node(value: node, side: .middle)
                let xL = Node(value: node-space.halfed, side: .right)
                let xR = Node(value: node+space.halfed, side: .left)
                let xP  = Node(value: node-space.step, side: .middle)
                let xLP = Node(value: node-space.step-space.halfed, side: .right)
                let xRP = Node(value: node-space.step+space.halfed, side: .left)
                guard let yL = Nf(t: tP, x: xL.value-c*time.step, xL: xLP, xM: xP, xR: xRP),
                      let yR = Nf(t: tP, x: xR.value-c*time.step, xL: xL, xM: x, xR: xR)
                else { return  }
                let z = u(xL.value, t)
                let v = u(xR.value, t)
                let n = u((xL.value+xR.value)/2, t)
                let y = 1.0/6.0*(z+4.0*n+v)
                detailed[t]?[x] = y
//                detailed[t]?[xL] = yL
//                detailed[t]?[xR] = yR
//                return
                guard (yR-y)*(y-yL) > 0 else {
                    detailed[t]?[xL] = y
                    detailed[t]?[xR] = y
                    return
                }
                let delta = yR-yL
                let sixth = 6.0*(y-0.5*(yL+yR))
                let deltaSix = delta*sixth
                let deltaSq = delta*delta
                detailed[t]?[xL] = deltaSix > deltaSq ? (3.0*y-2.0*yR) : yL
                detailed[t]?[xR] = deltaSix < -deltaSq ? (3.0*y-2.0*yL) : yR
            }
        }
    }
}
