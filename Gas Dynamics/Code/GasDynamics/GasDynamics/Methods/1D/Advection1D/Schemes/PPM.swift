//
//  PPM.swift
//  GasDynamics
//
//  Created by Арсений Токарев on 02.06.2022.
//

import Foundation

extension Advection1D {
    public final class PPM: Advection1D {
        public override func solve() {
            save(file: "advectionPPM", data: data)
        }
    }
}
