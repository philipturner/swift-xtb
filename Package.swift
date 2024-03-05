// swift-tools-version: 5.9
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
    .macOS(.v14)
  ],
  products: [
    .library(
      // Drop-in replacement for the JIT compiler that supports async copies.
      // TODO: Gate this module under an #if APPLE macro.
      name: "MetalCompiler",
      targets: ["MetalCompiler"]),
    .library(
      name: "Mechanosynthesis",
      targets: ["Mechanosynthesis"]),
    .library(
      name: "libxc",
      targets: ["libxc"]),
  ],
  dependencies: [
    // The OpenCL dependency is primarily for AMD and Nvidia GPUs.
    .package(url: "https://github.com/philipturner/swift-opencl", branch: "main"),
    .package(url: "https://github.com/philipturner/swift-numerics", branch: "Quaternions"),
  ],
  targets: [
    // Modules
    .target(
      name: "MetalCompiler",
      dependencies: []),
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
    
    // Tests
    .testTarget(
      name: "MetalCompilerTests",
      dependencies: [
        "MetalCompiler"
      ]),
    .testTarget(
      name: "MechanosynthesisTests",
      dependencies: [
        "Mechanosynthesis",
      ],
      cSettings: [
        CSetting.define("ACCELERATE_NEW_LAPACK")
      ]),
  ]
)
