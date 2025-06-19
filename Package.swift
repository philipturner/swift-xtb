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
        // WARNING: This binary appears to get copied into the '.build'
        // folder. It can cause troublesome, subtle problems where the dylib
        // gets recompiled, but the program doesn't register the change.
        //
        // In addition, it keeps referencing the dylib from Homebrew Cellar,
        // not the one pasted into the project's directory. Working on figuring
        // out why this goes wrong.
        //
        // I found a modification to the build workflow that solves both of
        // these problems simultaneously!
        .unsafeFlags(["-L\(FileManager.default.currentDirectoryPath)"]),
        .linkedLibrary("xtb.6")
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
