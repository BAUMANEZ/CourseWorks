//
//  PPM.swift
//  GasDynamics
//
//  Created by Арсений Токарев on 29.05.2022.
//

import Foundation

public class PPM1D: Algorithm1D {
    public func f(x: Double, t: Double) -> Double? {
        guard let mesh = solution[t],
              let cell = mesh.first(where: {
                  ($0.left?.x ?? x+1) <= x && ($0.right?.x ?? x-1) >= x
              }),
              let yL = cell.left?.y,
              let xi = xi(for: x, in: cell),
              let delta = delta(for: cell),
              let sixth = sixth(for: cell)
        else { return nil }
        return yL+xi*(delta+sixth*(1-x))
    }
    
    public override func solve() {
        super.solve()
        for t in time.nodes {
            guard var mesh = solution[t] else { continue }
            for (i, cell) in mesh.enumerated() {
                guard i > 0 && i < mesh.count-1 else { continue }
                let cellP = mesh[i+1]
                let cellM = mesh[i-1]
                var cell = cell
                cell.left?.y = cellM.right?.y ?? cellM.middle.y
                guard let y = cell.middle.y,
                      let yP = cellP.middle.y,
                      let delta = deltaM(for: cell, in: mesh),
                      let deltaP = deltaM(for: cellP, in: mesh)
                else { continue }
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
