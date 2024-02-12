// swift-tools-version: 5.8
// The swift-tools-version declares the minimum version of Swift required to build this package.

import PackageDescription

let package = Package(
  name: "DFT",
  products: [
    // Products define the executables and libraries a package produces, making them visible to other packages.
    .library(
      name: "DFT",
      targets: ["DFT"]),
  ],
  dependencies: [
    .package(url: "https://github.com/philipturner/swift-opencl", branch: "main"),
  ],
  targets: [
    // Targets are the basic building blocks of a package, defining a module or a test suite.
    // Targets can depend on other targets in this package and products from dependencies.
    .systemLibrary(
      name: "LibXC",
      pkgConfig: "libxc",
      providers: [
        .brew(["libxc"]),
        .apt(["libxc-dev"])
        // Not sure what to do on Windows yet.
      ]),
    .target(
      name: "DFT",
      dependencies: [
        .product(name: "OpenCL", package: "swift-opencl"),
        "LibXC",
      ]),
    .testTarget(
      name: "DFTTests",
      dependencies: ["DFT"]),
    .testTarget(
      name: "LibXCTests",
      dependencies: ["LibXC"]),
  ]
)
