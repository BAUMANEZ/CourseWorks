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
        public let steps : Int
        
        public final var step: Double {
            return (end-start)/Double(steps+1)
        }
        
        private var _nodes: [Double]?
        public final var nodes: [Double] {
            guard let cached = _nodes else {
                var values: [Double] = []
                values.reserveCapacity(steps)
                for x in stride(from: start, through: end, by: step) {
                    values.append(x)
                }
                _nodes = values
                return values
            }
            return cached
        }
        
        public init(start: Double, end: Double, steps: Int) {
            self.start = start
            self.end = end
            self.steps = steps
        }
        public init(start: Double, end: Double, step: Double) {
            self.start = start
            self.end = end
            self.steps = max(0, Int((end-start)/step) - 1)
        }
        
        public final func node(for i: Int) -> Double? {
            let node = Double(i)*step
            guard node >= start && node <= end else { return nil }
            return node
        }
        
        public static func == (lhs: Algorithm.Grid, rhs: Algorithm.Grid) -> Bool {
            return lhs.start == rhs.start && lhs.end == rhs.end && lhs.steps == rhs.steps
        }
    }
}

extension Algorithm1D.Grid {
    public static let empty = Algorithm.Grid(start: 0.0, end: 0.0, step: 0.0)
}
