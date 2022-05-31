//
//  PPM.swift
//  GasDynamics
//
//  Created by Арсений Токарев on 29.05.2022.
//

import Foundation

public final class PPM1D: Interpolation1D {
    public override func solve() {
        super.solve()
        for t in time.nodes {
            guard let mesh = solution[t] else { continue }
            var updated: Mesh = []
            for (i, cell) in mesh.enumerated() {
                var cell = cell
                cell.left?.y = updated.last?.right?.y
                guard i < mesh.count-1 else { continue }
                let cellP = mesh[i+1]
                guard let y = cell.middle.y,
                      let yP = cellP.middle.y,
                      let delta = deltaM(for: cell, in: mesh),
                      let deltaP = deltaM(for: cellP, in: mesh)
                else { continue }
                cell.right?.y = 0.5*(y+yP)-1/6*(deltaP-delta)
                updated.append(cell)
            }
            adjust(mesh: &updated)
            solution[t] = updated
        }
        save(file: "advectionPPM", data: interpolationData)
    }
}

//MARK: - Helpers
extension PPM1D {
    public final func deltaM(for cell: Cell, in mesh: Mesh) -> Double? {
        guard let y = cell.middle.y,
              let i = mesh.firstIndex(of: cell)
        else { return nil }
        guard i > 0 else {
            guard mesh.indices.contains(i+1),
                  let yP = mesh[i+1].middle.y
            else { return nil }
            guard yP-y > 0 else { return .zero }
            let delta = 0.5*(yP+y)
            let min1 = 2*fabs(yP-y)
            return delta.sign*min(fabs(delta), min1)
        }
        guard i < mesh.count else {
            guard mesh.indices.contains(i-1),
                  let yM = mesh[i-1].middle.y
            else { return nil }
            let delta = 0.5*(y+yM)
            let min2 = 2*fabs(y-yM)
            return delta.sign*min(fabs(delta), min2)
        }
        guard mesh.indices.contains(i-1),
              let yM = mesh[i-1].middle.y,
              mesh.indices.contains(i+1),
              let yP = mesh[i+1].middle.y
        else { return nil }
        let right = yP-y
        let left = y-yM
        guard right*left > 0 else { return 0.0 }
        let delta = 0.5*(yP+yM)
        let min1 = 2*fabs(right)
        let min2 = 2*fabs(left)
        return delta.sign*min(fabs(delta), min1, min2)
    }
}

extension Double {
    public var sign: Double {
        self >= 0 ? 1.0 : -1.0
    }
}
