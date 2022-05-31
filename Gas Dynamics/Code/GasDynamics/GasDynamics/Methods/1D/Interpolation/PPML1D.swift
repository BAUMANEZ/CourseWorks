//
//  PPML.swift
//  GasDynamics
//
//  Created by Арсений Токарев on 29.05.2022.
//

import Foundation

public final class PPML1D: Interpolation1D {
    public let c: () -> Double
    
    public required init(initial: Algorithm1D) {
        if let algorithm = initial as? Advection1D {
            self.c = { return algorithm.c }
        } else if let algorithm = initial as? Burgers1D {
            self.c = { return .zero }
        } else {
            self.c = { return .zero }
        }
        super.init(initial: initial)
    }
    
    public required init(a: Double, b: Double, h: Double, tau: Double, deadline: Double) {
        fatalError("init(a:b:h:tau:deadline:) has not been implemented")
    }
    
    public required init(tau: Double = 1.0, deadline: Double) {
        fatalError("init(tau:deadline:) has not been implemented")
    }
    
    public override func solve() {
        guard time.nodes.count > 1 else { return }
        for t in time.nodes[1...] {
            guard let mesh = solution[t],
                  let prevMesh = solution[t-time.step]
            else { continue }
            var updated: Mesh = []
            let c = c()
            if c > 0 {
                for (i, cellM) in prevMesh.enumerated() {
                    var cell = mesh[i]
                    cell.left?.y = updated.last?.right?.y
                    guard let delta = delta(for: cellM),
                          let sixth = sixth(for: cellM),
                          let yL = cellM.left?.y ?? cellM.middle.y
                    else { continue }
                    let xi = 1-c*time.step/space.step
                    cell.right?.y = yL+xi*(delta+sixth*(1-xi))
                    updated.append(cell)
                }
            }
            adjust(mesh: &updated)
            solution[t] = updated
        }
        save(file: "advectionPPML", data: interpolationData)
    }
}
