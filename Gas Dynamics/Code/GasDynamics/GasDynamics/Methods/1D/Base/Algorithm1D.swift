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
    
    public var plotStep: Int {
        return 40
    }
    public init(a: Double, b: Double, h: Double, tau: Double, deadline: Double) {
        self.space = Grid(start: a, end: b, step: h)
        super.init(tau: tau, deadline: deadline)
    }
    
    public final func data(for solutions: [Time: Mesh]) -> Data? {
        let timestamps = solutions.keys.sorted(by: { $0 < $1 })
        let json = solutions.reduce(into: [String: Any]()) { json, solution in
            guard let index = timestamps.firstIndex(of: solution.key),
                  index%plotStep == 0
            else { return () }
            let mesh = solution.value.reduce(into: [String: String]()) { mesh, pair in
                let x = pair.key
                let y = pair.value
                mesh[String(x)] = String(y)
            }
            json[String(solution.key)] = mesh
        }
        return try? JSONSerialization.data(withJSONObject: json, options: .sortedKeys)
    }
    
    public func f(x: Double, t: Double) -> Double? {
        return nil
    }
}
