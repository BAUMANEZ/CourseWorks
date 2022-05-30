//
//  PPM.swift
//  GasDynamics
//
//  Created by Арсений Токарев on 29.05.2022.
//

import Foundation

public class PPM1D: Algorithm1D {
    public override func solve() {
        super.solve()
        for t in time.nodes {
            guard var mesh = solution[t] else { continue }
            for (segment, value) in mesh {
                guard let index = segments.firstIndex(of: segment),
                      index > 0 || index < segments.count-1,
                      let segmentP = self.segment(for: index+1),
                      let yP = self.value(for: segmentP, in: mesh),
                      let deltaP = deltaM(for: segmentP, in: mesh),
                      let delta = deltaM(for: segment, in: mesh)
                else { continue }
                var value = value
                value.left = self.value(for: self.segment(for: index-1), in: mesh)?.right
                value.right = 0.5*(value.middle+yP.middle)-1/6*(deltaP-delta)
                mesh[segment] = value
            }
            solution[t] = mesh
        }
        save(file: "advectionPPM")
    }
}

//MARK: - Helpers
extension PPM1D {
    public func deltaM(for segment: Boundary, in mesh: Mesh) -> Double? {
        guard let values = mesh[segment],
              let index = segments.firstIndex(of: segment),
              let valuesP = value(for: self.segment(for: index+1), in: mesh),
              let valuesM = value(for: self.segment(for: index-1), in: mesh)
        else{ return nil }
        let right = valuesP.middle-values.middle
        let left = values.middle-valuesM.middle
        guard right*left > 0 else { return 0.0 }
        let delta = 0.5*(valuesP.middle+valuesM.middle)
        let min1 = 2*fabs(right)
        let min2 = 2*fabs(left)
        return delta.sign*min(delta, min1, min2)
    }
}

extension Double {
    public var sign: Double {
        self >= 0 ? 1.0 : -1.0
    }
}
