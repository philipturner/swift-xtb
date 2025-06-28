// swift-tools-version: 6.1

import PackageDescription

var linkerSettings: [LinkerSetting] = []

// Best attempt at scoping this peculiar search path to macOS for now. MM4
// handles it differently, with an environment variable for the path.
#if os(macOS)
import class Foundation.FileManager

linkerSettings += [
  .unsafeFlags(["-L\(FileManager.default.currentDirectoryPath)"]),
  .linkedLibrary("xtb_accelerate") // change to 'xtb' before supporting Windows
]
#endif

let package = Package(
  name: "swift-xtb",
  dependencies: [],
  targets: [
    .target(
      name: "C_xTB",
      dependencies: [],
      linkerSettings: linkerSettings),
    .target(
      name: "xTB",
      dependencies: ["C_xTB"]),
    .executableTarget(
      name: "Workspace",
      dependencies: ["xTB"]),
  ]
)
