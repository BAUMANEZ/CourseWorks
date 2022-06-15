//
//  main.swift
//  GasDynamics
//
//  Created by Арсений Токарев on 29.05.2022.
//

import Foundation

let c: Double = 1.0
let h: Double = 1.0/2
let profile: Advection1D.Profile = .tooth

//let ppm = Advection1D.PPM(c: c, h: h, profile: profile)
let ppml = Advection1D.PPML(c: c, h: h, profile: profile)

//for h in [1.0, 0.5, 0.25, 0.125, 0.0625] {
//    print(Advection1D.PPML(c: c, h: h, profile: .cos).normL2)
//}
