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
                var yL: Double?
                if let _yL = detailed[t]?[Node(value: xL.value, side: .left)] {
                    yL = _yL
                    detailed[t]?[xL] = _yL
                }
                let xM = Node(value: node-space.step, side: .middle)
                let xP = Node(value: node+space.step, side: .middle)
                let xPP = Node(value: node+2.0*space.step, side: .middle)
                guard let y   = detailed[t]?[x],
                      let yM  = detailed[t]?[xM],
                      let yP  = detailed[t]?[xP],
                      let yPP = detailed[t]?[xPP]
                else { return () }
                let yR = 0.5*(y+yP)-(1.0/6.0)*(deltaM(yL: y, y: yP, yR: yPP)-deltaM(yL: yM, y: y, yR: yP))
                detailed[t]?[xR] = yR
                if let yL = yL {
                    guard (yR-y)*(y-yL) > 0 else {
                        detailed[t]?[xL] = y
                        detailed[t]?[xR] = y
                        return
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
