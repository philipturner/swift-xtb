// swift-tools-version: 5.9
// The swift-tools-version declares the minimum version of Swift required to build this package.

import PackageDescription

let package = Package(
  name: "swift-xtb",
  products: [
    .library(
      name: "xTB",
      targets: ["xTB"]),
  ],
  dependencies: [],
  targets: [
    .target(
      name: "xTB",
      dependencies: []),
  ]
)
