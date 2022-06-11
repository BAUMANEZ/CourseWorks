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
                let xMM = Node(value: node-2.0*space.step, side: .middle)
                let xM  = Node(value: node-space.step, side: .middle)
                let xP  = Node(value: node+space.step, side: .middle)
                let xPP = Node(value: node+2.0*space.step, side: .middle)
                guard let yMM = detailed[t]?[xMM],
                      let yM  = detailed[t]?[xM],
                      let y   = detailed[t]?[x],
                      let yP  = detailed[t]?[xP],
                      let yPP = detailed[t]?[xPP]
                else { return () }
                let xL = Node(value: node-space.halfed, side: .right)
                let xR = Node(value: node+space.halfed, side: .left)
                let yL = 0.5*(yM+y)-(1.0/6.0)*(deltaM(yL: yM, y: y, yR: yP)-deltaM(yL: yMM, y: yM, yR: y))
                let yR = 0.5*(y+yP)-(1.0/6.0)*(deltaM(yL: y, y: yP, yR: yPP)-deltaM(yL: yM, y: y, yR: yP))
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
        
        private func deltaM(yL: Double, y: Double, yR: Double) -> Double {
            let deltaM = 0.5*(yL+yR)
            guard ((yR-y)*(y-yL)).sign > 0 else { return .zero }
            let min1 = 2.0*(yR-y)
            let min2 = 2.0*(y-yL)
            return deltaM.sign*min(Swift.abs(deltaM), Swift.abs(min1), Swift.abs(min2))
        }
    }
}

public extension Double {
    var sign: Double {
        return self >= .zero ? 1.0 : -1.0
    }
}
