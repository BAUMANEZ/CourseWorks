//
//  main.swift
//  GasDynamics
//
//  Created by Арсений Токарев on 29.05.2022.
//

import Foundation

let profile: Advection1D.Profile = .tooth
let advection = Advection1D(c: 1.0, h: 1.0, profile: profile)
let ppm = Advection1D.PPM(c: 1.0, h: 1.0, profile: profile)
let ppml = Advection1D.PPML(c: 1.0, h: 1.0, profile: profile)
