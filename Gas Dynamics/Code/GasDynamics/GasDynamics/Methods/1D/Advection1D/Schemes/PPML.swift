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
                if let yL = detailed[t]?[Node(value: xL.value, side: .left)] {
                    detailed[t]?[xL] = yL
                }
                guard let yLP = detailed[tP]?[xL],
                      let yP  = detailed[tP]?[x],
                      let yRP = detailed[tP]?[xR]
                else { return () }
                let xi = 1.0-gamma
                let delta = yRP-yLP
                let sixth = 6.0*(yP-0.5*(yLP+yRP))
                detailed[t]?[xR] = yLP+xi*(delta+sixth*(1.0-xi))
            }
        }
    }
}
