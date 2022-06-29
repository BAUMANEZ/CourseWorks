//
//  main.swift
//  GasDynamics
//
//  Created by Арсений Токарев on 29.05.2022.
//

import Foundation

let c: Double = 1.0
let h: Double = 1.0/4.0
let sigma: Double = 1.0/4.0
let profile: Advection1D.Profile = .cos

let method = Advection1D.PPML(c: c, h: h, sigma: sigma, profile: profile)
print(method.normC)
print(method.normL)
print(method.normL2)
