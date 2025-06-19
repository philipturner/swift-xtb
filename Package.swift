// swift-tools-version: 5.9
// The swift-tools-version declares the minimum version of Swift required to build this package.

import PackageDescription
import class Foundation.FileManager

let package = Package(
  name: "swift-xtb",
  dependencies: [],
  targets: [
    .target(
      name: "C_xTB",
      dependencies: [],
      linkerSettings: [
        .unsafeFlags(["-L\(FileManager.default.currentDirectoryPath)"]),
        .linkedLibrary("xtb_accelerate")
      ]),
    .target(
      name: "xTB",
      dependencies: [
        "C_xTB"
      ]),
    .executableTarget(
      name: "Workspace",
      dependencies: [
        "xTB"
      ]),
  ]
)
