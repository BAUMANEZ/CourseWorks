//
//  main.swift
//  GasDynamics
//
//  Created by Арсений Токарев on 29.05.2022.
//

import Foundation

let c: Double = 1.0
let h: Double = 1.0
let profile: Advection1D.Profile = .cos

//let ppm = Advection1D.PPM(c: c, h: h, profile: profile)
//let ppml = Advection1D.PPML(c: c, h: h, profile: profile)

print("PPM")
for h in [1.0, 0.5, 0.25, 0.125, 0.0625] {
    print(Advection1D.PPM(c: c, h: h, profile: profile).normL2)
}
print("----")
print("PPML")
for h in [1.0, 0.5, 0.25, 0.125, 0.0625] {
    print(Advection1D.PPML(c: c, h: h, profile: profile).normL2)
}
