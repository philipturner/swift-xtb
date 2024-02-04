// swift-tools-version: 5.9
// The swift-tools-version declares the minimum version of Swift required to build this package.

import PackageDescription

let package = Package(
    name: "DFT",
    platforms: [
      // Specifying a macOS/iOS version doesn't affect deployment to non-Apple platforms.
      .macOS(.v13),
      .iOS(.v15),
    ],
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
        .target(
            name: "DFT",
            dependencies: [
              .product(name: "OpenCL", package: "swift-opencl"),
            ]),
        .testTarget(
            name: "DFTTests",
            dependencies: ["DFT"]),
    ]
)
