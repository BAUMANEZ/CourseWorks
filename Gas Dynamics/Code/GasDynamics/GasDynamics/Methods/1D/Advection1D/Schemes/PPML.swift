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
        
        public override func solve(for j: Int) {
            super.solve(for: j)
            guard detailed[j] != nil, detailed[j-1] != nil else { return }
            space.nodes().forEach { node in
                let x  = Node(value: node, side: .middle)
                let xL = Node(value: node-space.halfed, side: .right)
                let xR = Node(value: node+space.halfed, side: .left)
                let xP  = Node(value: node-space.step, side: .middle)
                let xLP = Node(value: node-space.step-space.halfed, side: .right)
                let xRP = Node(value: node-space.halfed, side: .left)
                guard let yL = Nf(j: j-1, x: xL.value-c*time.step, xL: xLP, xM: xP, xR: xRP),
                      let yR = Nf(j: j-1, x: xR.value-c*time.step, xL: xL, xM: x, xR: xR),
                      let y  = detailed[j]?[x]
                else { return  }
                guard (yR-y)*(y-yL) > 0 else {
                    detailed[j]?[xL] = y
                    detailed[j]?[xR] = y
                    return
                }
                let delta = yR-yL
                let sixth = 6.0*(y-0.5*(yL+yR))
                let deltaSix = delta*sixth
                let deltaSq = delta*delta
                detailed[j]?[xL] = deltaSix > deltaSq ? (3.0*y-2.0*yR) : yL
                detailed[j]?[xR] = deltaSix < -deltaSq ? (3.0*y-2.0*yL) : yR
            }
        }
    }
}
