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
            guard let mesh = detailed[t] else { return }
            detailed[t] = space.nodes.reduce(into: DetailedMesh()) { detailed, node in
                let x = Node(value: node, side: .middle)
                guard let y = mesh[x] else { return () }
                detailed[x] = y
                let nodeL = Node(value: node-space.halfed, side: .left)
                if let yL = detailed[nodeL] {
                    detailed[Node(value: nodeL.value, side: .right)] = yL
                }
                let xM = Node(value: node-space.step, side: .middle)
                let xP = Node(value: node+space.step, side: .middle)
                let xPP = Node(value: node+2.0*space.step, side: .middle)
                guard let yM = mesh[xM],
                      let yP = mesh[xP],
                      let yPP = mesh[xPP]
                else { return () }
                let xR = Node(value: node+space.halfed, side: .left)
                detailed[xR] = 0.5*(y+yP)-1/6*(deltaM(yL: y, y: yP, yR: yPP)-deltaM(yL: yM, y: y, yR: yP))
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
