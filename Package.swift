// swift-tools-version: 5.9
// The swift-tools-version declares the minimum version of Swift required to build this package.

import PackageDescription

let package = Package(
  name: "Mechanosynthesis",
  products: [
    .library(
      name: "Mechanosynthesis",
      targets: ["Mechanosynthesis"]),
  ],
  dependencies: [
    .package(url: "https://github.com/philipturner/swift-numerics", branch: "Quaternions"),
  ],
  targets: [
    .target(
      name: "Mechanosynthesis",
      dependencies: [
        .product(name: "Numerics", package: "swift-numerics"),
      ]),
    .testTarget(
      name: "MechanosynthesisTests",
      dependencies: ["Mechanosynthesis"]),
  ]
)
