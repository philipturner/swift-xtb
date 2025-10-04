// swift-tools-version: 6.1

import PackageDescription

let package = Package(
  name: "swift-xtb",
  dependencies: [],
  targets: [
    .target(
      name: "xTB",
      dependencies: []),
    .executableTarget(
      name: "Workspace",
      dependencies: ["xTB"]),
  ]
)
