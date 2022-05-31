//
//  main.swift
//  GasDynamics
//
//  Created by Арсений Токарев on 29.05.2022.
//

import Foundation

let advection = Advection1D(profile: .leftTriangle)
let ppm = PPM1D(initial: advection)
let ppml = PPML1D(initial: advection)
//print(ppm.normC)
