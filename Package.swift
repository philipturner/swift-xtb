// swift-tools-version: 5.9
// The swift-tools-version declares the minimum version of Swift required to build this package.

import PackageDescription

let package = Package(
  name: "Mechanosynthesis",
  products: [
    .library(
      name: "LinearAlgebra",
      targets: ["LinearAlgebra"]),
    .library(
      name: "Meshing",
      targets: ["Meshing"]),
  ],
  dependencies: [
    .package(url: "https://github.com/philipturner/swift-numerics", branch: "Quaternions"),
  ],
  targets: [
    .target(
      name: "LinearAlgebra",
      dependencies: []),
    .target(
      name: "Meshing",
      dependencies: [
        .product(name: "Numerics", package: "swift-numerics"),
      ]),
    .target(
      name: "xTB",
      dependencies: []),
    
    .executableTarget(
      name: "Workspace",
      dependencies: [
        "LinearAlgebra", 
        "Meshing",
        .product(name: "Numerics", package: "swift-numerics"),
        "xTB",
      ]),
    
    .testTarget(
      name: "AlgorithmTests",
      dependencies: [
        "LinearAlgebra",
        .product(name: "Numerics", package: "swift-numerics"),
      ]),
    .testTarget(
      name: "MeshingTests",
      dependencies: [
        "Meshing"
      ]),
  ]
)
