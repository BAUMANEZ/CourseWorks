//
//  Grid.swift
//  GasDynamics
//
//  Created by Арсений Токарев on 02.06.2022.
//

import Foundation

extension Algorithm {
    public final class Grid: Equatable {
        //MARK: Stepper
        /// - parameter:  order of node (Int)
        /// - result: value of node (Double)
        public typealias StepProvider = (Int) -> Double
        
        public let start : Double
        public let end   : Double
        public let step  : Double
        public let steps : Int
        
        public final func range(starting: Int? = nil, ending: Int? = nil) -> ClosedRange<Int> {
            let starting = starting ?? 0
            var ending = ending ?? steps
            if ending < 0 { ending = steps-ending }
            guard ending >= starting else { return Int.max...Int.max }
            return starting...ending
        }
        
        public final func nodes(starting: Int? = nil, ending: Int? = nil) -> [Double] {
            return range(starting: starting, ending: ending).compactMap{ self.node(for: $0) }
        }
        
        public final var halfed: Double {
            return step/2
        }
        
        public final var json: [String: String] {
            return [
                "a": String(start),
                "b": String(end),
                "h": String(step),
                "n": String(steps),
            ]
        }
        
        
        public init(start: Double, end: Double, step: Double) {
            self.start = start
            self.end = end
            self.step = step
            self.steps = max(0, Int((end-start)/step) - 1)
        }
        
        public final func node(for i: Int) -> Double {
            return start+Double(i)*step
        }
        
        public static func == (lhs: Algorithm.Grid, rhs: Algorithm.Grid) -> Bool {
            return lhs.start == rhs.start && lhs.end == rhs.end && lhs.steps == rhs.steps
        }
    }
}

extension Algorithm1D.Grid {
    public static let empty = Algorithm.Grid(start: 0.0, end: 0.0, step: 0.0)
}
