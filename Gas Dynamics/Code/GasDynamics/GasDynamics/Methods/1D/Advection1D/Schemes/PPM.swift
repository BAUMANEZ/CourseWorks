//
//  PPM.swift
//  GasDynamics
//
//  Created by Арсений Токарев on 02.06.2022.
//

import Foundation

extension Advection1D {
    public final class PPM: InterpolatedAdvection1D {
        public override var identifier: String? {
            return "advectionPPM"
        }
        
        public override func solve(for t: Time) {
            super.solve(for: t)
            guard detailed[t] != nil else { return }
            space.nodes.forEach { node in
                let x  = Node(value: node, side: .middle)
                let xL = Node(value: node-space.halfed, side: .right)
                let xR = Node(value: node+space.halfed, side: .left)
                if let yL = detailed[t]?[Node(value: xL.value, side: .left)] {
                    detailed[t]?[xL] = yL
                }
                let xM = Node(value: node-space.step, side: .middle)
                let xP = Node(value: node+space.step, side: .middle)
                let xPP = Node(value: node+2.0*space.step, side: .middle)
                guard let y   = detailed[t]?[x],
                      let yM  = detailed[t]?[xM],
                      let yP  = detailed[t]?[xP],
                      let yPP = detailed[t]?[xPP]
                else { return () }
                detailed[t]?[xR] = 0.5*(y+yP)-(1.0/6.0)*(deltaM(yL: y, y: yP, yR: yPP)-deltaM(yL: yM, y: y, yR: yP))
            }
        }
        
        private func deltaM(yL: Double, y: Double, yR: Double) -> Double {
            let deltaM = 0.5*(yL+yR)
            guard (yR-y)*(y-yL) > .zero else { return .zero }
            let min1 = 2.0*(yR-y)
            let min2 = 2.0*(y-yL)
            return deltaM.sign*min(fabs(deltaM), fabs(min1), fabs(min2))
        }
    }
}

public extension Double {
    var sign: Double {
        return self >= .zero ? 1.0 : -1.0
    }
}
