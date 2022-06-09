//
//  ArticleProfiles.swift
//  GasDynamics
//
//  Created by Арсений Токарев on 02.06.2022.
//

import Foundation

extension Advection1D {
    public enum Profile {
        case leftTriangle
        case rightTriangle
        case rectangle
        case tooth
        case M
        case cos
        
        public var l1 : Double { return  10.0 }
        public var l11: Double { return  50.0/3.0 }
        public var l12: Double { return  20.0 }
        public var l22: Double { return  70.0/3.0 }
        public var l2 : Double { return  30.0 }
        public var L  : Double { return  400.0 }
        public var T  : Double { return  400.0 }
        
        public func f(x: Double) -> Double {
            guard x >= l1 && x <= l2 else { return .zero }
            switch self {
            case .leftTriangle:
                return (x-l1)/(l2-l1)
            case .rightTriangle:
                return (l2-x)/(l2-l1)
            case .rectangle:
                return 1.0
            case .cos:
                return 0.5-0.5*Foundation.cos(2.0*Double.pi/(l2-l1) * (x-l1) )
            case .tooth:
                switch x {
                case l1 ..< l11:
                    return -2.0*(x-l1)/(3.0*(l11-l1)) + 1.0
                case l11 ... l22:
                    return 1.0/3.0
                default:
                    return 2.0*(x-l2)/(3.0*(l2-l22)) + 1.0
                }
            case .M:
                switch x {
                case l1 ..< l12:
                    return -2.0*(x-l1)/(3.0*(l12-l1)) + 1.0
                default:
                    return 2.0*(x-l2)/(3.0*(l2-l12)) + 1.0
                }
            }
        }
        
        public var description: String {
            switch self {
            case .rightTriangle: return "rightTriangle"
            case .leftTriangle : return "leftTriangle"
            case .rectangle    : return "rectangle"
            case .tooth        : return "tooth"
            case .cos          : return "cos"
            case .M            : return "m"
            }
        }
    }
}
