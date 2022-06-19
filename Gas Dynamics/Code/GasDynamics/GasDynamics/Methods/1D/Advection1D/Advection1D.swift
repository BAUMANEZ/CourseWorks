//
//  Advection.swift
//  GasDynamics
//
//  Created by Арсений Токарев on 31.05.2022.
//

import Foundation

public class Advection1D: Algorithm1D {
    public let c : Double
    public let u : (Double, Time) -> Double
    public let u0: (Double) -> Double
    public let profile: String?
    
    public var gamma: Double {
        return c*time.step/space.step
    }
    
    public init(a : Double,
                b : Double,
                T : Double,
                c : Double,
                h : Double,
                u : @escaping (Double, Time) -> Double,
                u0: @escaping (Double) -> Double,
                profile: String?
    ) {
        self.c = c
        self.u = u
        self.u0 = u0
        self.profile = profile
        super.init(a: a, b: b, h: h, tau: 0.5*h/c, deadline: T)
    }
    
    public convenience init(c: Double = 1.0, h: Double = 1.0, profile: Profile) {
        let u = { x, t in return profile.f(x: x-c*t) }
        self.init(a: profile.l1, b: profile.L, T: profile.T, c: c, h: h, u: u, u0: profile.f, profile: profile.description)
    }
    
    //MARK: Drift along characteristics
    public final func drift(for t: Time, in x: Double) -> Double {
        return u0(x-c*t)
    }
    public override func f(x: Double, t: Double) -> Double? {
        return drift(for: t, in: x)
    }
    
    public override func solve() {
        let solutions = time.nodes(starting: 1).reduce(into: [Time: Mesh]()) { solutions, t in
            solutions[t] = space.nodes().reduce(into: Mesh()) { mesh, x in
                mesh[x] = f(x: x, t: t)
            }
        }
        save(file: "advection", data: data(for: solutions))
    }
}
