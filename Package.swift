// swift-tools-version: 6.1

import PackageDescription

let package = Package(
  name: "swift-xtb",
  dependencies: [
    .package(url: "https://github.com/apple/swift-docc-plugin", branch: "main"),
  ],
  targets: [
    .target(
      name: "xTB",
      dependencies: []),
    .testTarget(
      name: "xTB_Tests",
      dependencies: ["xTB"]),
  ]
)
