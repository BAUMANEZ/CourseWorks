//
//  Base.swift
//  GasDynamics
//
//  Created by Арсений Токарев on 31.05.2022.
//

import Foundation

public class Algorithm {
    public let time: Grid
    
    public var data: Data? {
        return nil
    }
    
    public init(tau: Double = 1.0, deadline: Double) {
        self.time = Grid(start: 0, end: deadline, step: tau)
        solve()
    }
    
    public func solve() {}
    
    public final func save(file name: String, data: Data?, meta: Data? = nil) {
        let manager = FileManager.default
        let folder = manager.homeDirectoryForCurrentUser.appendingPathComponent("Desktop", isDirectory: true).appendingPathComponent("CourseWorks", isDirectory: true).appendingPathComponent("Gas Dynamics", isDirectory: true).appendingPathComponent("Code", isDirectory: true).appendingPathComponent("data", isDirectory: true)
        let file = folder.appendingPathComponent(name).appendingPathExtension("json")
        let extra = folder.appendingPathComponent("metadata").appendingPathExtension("json")
        guard let _ = try? manager.createDirectory(at: folder, withIntermediateDirectories: true) else { return }
        manager.createFile(atPath: file.path, contents: data)
        manager.createFile(atPath: extra.path, contents: meta)
    }
}
