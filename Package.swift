// swift-tools-version: 5.8
// The swift-tools-version declares the minimum version of Swift required to build this package.

import PackageDescription
import class Foundation.ProcessInfo

var linkerSettings: [LinkerSetting] = []

if let path = ProcessInfo.processInfo.environment["XC_LIBRARY_PATH"] {
  linkerSettings = [
    .unsafeFlags(["-L\(path)"]),
    .linkedLibrary("xc"),
  ]
}

let package = Package(
  name: "mechanosynthesis",
  platforms: [
    .macOS(.v13)
  ],
  products: [
    // Products define the executables and libraries a package produces, making them visible to other packages.
    .library(
      name: "Mechanosynthesis",
      targets: ["Mechanosynthesis"]),
    .library(
      name: "libxc",
      targets: ["libxc"]),
  ],
  dependencies: [
    .package(url: "https://github.com/philipturner/swift-opencl", branch: "main"),
    .package(url: "https://github.com/philipturner/swift-numerics", branch: "Quaternions"),
  ],
  targets: [
    // Targets are the basic building blocks of a package, defining a module or a test suite.
    // Targets can depend on other targets in this package and products from dependencies.
    .target(
      name: "Mechanosynthesis",
      dependencies: [
        "libxc",
        .product(name: "Numerics", package: "swift-numerics"),
        .product(name: "OpenCL", package: "swift-opencl"),
      ]),
    .target(
      name: "libxc",
      dependencies: [],
    linkerSettings: linkerSettings),
    .testTarget(
      name: "MechanosynthesisTests",
      dependencies: [
        "Mechanosynthesis",
      ]),
  ]
)
