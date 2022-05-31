//
//  Interpolation1D.swift
//  GasDynamics
//
//  Created by Арсений Токарев on 31.05.2022.
//

import Foundation

public class Interpolation1D: Algorithm1D {
    public required init(initial: Algorithm1D) {
        super.init(a: initial.space.start,
                   b: initial.space.end,
                   h: initial.space.step,
                   tau: initial.time.step,
                   deadline: initial.time.end)
        self.solution = initial.solution
        solve()
    }
    
    public required init(a: Double, b: Double, h: Double, tau: Double, deadline: Double) {
        super.init(a: a, b: b, h: h, tau: tau, deadline: deadline)
    }
    
    public required init(tau: Double = 1.0, deadline: Double) {
        fatalError("init(tau:deadline:) has not been implemented")
    }
}
