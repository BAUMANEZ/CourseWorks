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
            for (i, cell) in mesh.enumerated() {
                guard let y = cell.middle.y,
                      let index = mesh.firstIndex(of: cell),
                      index > 0 && index < mesh.count-1,
                      mesh.indices.contains(index+1),
                      let yP = mesh[index+1].middle.y,
                      let delta = deltaM(for: cell, in: mesh),
                      let deltaP = deltaM(for: mesh[index+1], in: mesh)
                else { continue }
                let cellP = mesh[index+1]
                var cell = cell
                cell.left = cellP.right
                cell.right?.y = 0.5*(y+yP)-1/6*(deltaP-delta)
                mesh[i] = cell
            }
            solution[t] = mesh
        }
        save(file: "advectionPPM")
    }
}

//MARK: - Helpers
extension PPM1D {
    public func deltaM(for cell: Cell, in mesh: Mesh) -> Double? {
        guard let y = cell.middle.y,
              let index = mesh.firstIndex(of: cell),
              mesh.indices.contains(index-1),
              let yM = mesh[index-1].middle.y,
              mesh.indices.contains(index+1),
              let yP = mesh[index+1].middle.y
        else { return nil }
        let right = yP-y
        let left = y-yM
        guard right*left > 0 else { return 0.0 }
        let delta = 0.5*(yP+yM)
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
