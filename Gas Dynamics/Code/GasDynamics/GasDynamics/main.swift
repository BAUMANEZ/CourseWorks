//
//  main.swift
//  GasDynamics
//
//  Created by Арсений Токарев on 29.05.2022.
//

import Foundation

let ppm = PPM1D(initial: Advection1D(profile: .cos))
print(ppm.f(x: ppm.space.start, t: 0))
