//
//  Algorithm.swift
//  GasDynamics
//
//  Created by Арсений Токарев on 29.05.2022.
//

import Foundation

public class Algorithm1D: Algorithm {
    public typealias Time = Double
    public typealias Mesh = Dictionary<Double, Double>
    
    public let space: Grid
    public var solutions: [Time: Mesh] = [:]
    
    public var plotStep: Int {
        return 35
    }
    
    public override var data: Data? {
        let json = solutions.reduce(into: [String: Any]()) { json, solution in
            guard let step = time.nodes.firstIndex(of: solution.key), step%plotStep == 0 else { return () }
            let mesh = solution.value.reduce(into: [String: String]()) { mesh, pair in
                let x = pair.key
                let y = pair.value
                mesh[String(x)] = String(y)
            }
            json[String(solution.key)] = mesh
        }
        return try? JSONSerialization.data(withJSONObject: json, options: .sortedKeys)
    }
    
    public init(a: Double, b: Double, h: Double, tau: Double, deadline: Double) {
        self.space = Grid(start: a, end: b, step: h)
        super.init(tau: tau, deadline: deadline)
    }
    
    public final func f(x: Double, t: Double) -> Double? {
        return nil
    }
}

//MARK: - Helpers
extension Algorithm1D {
}
