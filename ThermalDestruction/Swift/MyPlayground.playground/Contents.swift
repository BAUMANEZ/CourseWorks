import UIKit

func getDocumentsDirectory() -> URL {
    let paths = FileManager.default.homeDirectoryForCurrentUser
    return paths[0]
}

let str = "Super long string here"
let filename =  getDocumentsDirectory().appendingPathComponent("output.txt")

do {
    try str.write(to: filename, atomically: true, encoding: String.Encoding.utf8)
} catch {
    // failed to write file – bad permissions, bad filename, missing permissions, or more likely it can't be converted to the encoding
}
