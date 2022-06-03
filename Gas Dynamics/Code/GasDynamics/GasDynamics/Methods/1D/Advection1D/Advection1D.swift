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
    
    public var gamma: Double {
        return space.step*c/time.step
    }
    
    public init(a : Double,
                b : Double,
                T : Double,
                c : Double,
                h : Double,
                u : @escaping (Double, Time) -> Double,
                u0: @escaping (Double) -> Double
    ) {
        self.c = c
        self.u = u
        self.u0 = u0
        super.init(a: a, b: b, h: 1.0, tau: c/h, deadline: T)
    }
    
    public convenience init(c: Double = 1.0, h: Double = 1.0, profile: Profile) {
        let u = { x, t in return profile.f(x: x+c*t) }
        self.init(a: profile.l1, b: profile.L, T: profile.T, c: c, h: h, u: u, u0: profile.f)
    }
    
    public override func solve() {
        solutions[time.start] = space.nodes.reduce(into: Mesh()) { mesh, x in
            mesh[x] = self.u0(x)
        }
        if time.steps > 1 {
            let gamma = gamma
            for t in time.nodes[1...] {
                guard let meshP = solutions[t-time.step] else { continue }
                var mesh = Mesh()
                for x in space.nodes {
                    guard let yP = meshP[x], let yPM = meshP[x-space.step] else {
                        mesh[x] = 0.0
                        continue
                    }
                    mesh[x] = (1-gamma)*yP+gamma*yPM
                }
                solutions[t] = mesh
            }
        }
        save(file: "advection", data: data)
    }
}
