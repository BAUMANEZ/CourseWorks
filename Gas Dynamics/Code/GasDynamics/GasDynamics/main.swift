//
//  main.swift
//  GasDynamics
//
//  Created by Арсений Токарев on 29.05.2022.
//

import Foundation

let c = 1.0
let h = 1.0

let profile: Advection1D.Profile = .cos
let advection = Advection1D(c: c, h: h, profile: profile)
let ppm = Advection1D.PPM(c: c, h: h, profile: profile)
let ppml = Advection1D.PPML(c: c, h: h, profile: profile)
print(ppm.normC)
print(ppml.normC)
