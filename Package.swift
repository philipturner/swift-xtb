// swift-tools-version: 5.9
// The swift-tools-version declares the minimum version of Swift required to build this package.

import PackageDescription
import class Foundation.ProcessInfo

var targets: [Target] = []

// MARK: - Apple Platforms

var platforms: [SupportedPlatform] = []
#if canImport(Darwin)
platforms.append(.macOS(.v13))
#endif

var cSettings: [CSetting] = []
var swiftSettings: [SwiftSetting] = []

// MARK: - Package Dependencies

var packageDependencies: [Package.Dependency] = []
packageDependencies.append(
  .package(
    url: "https://github.com/philipturner/swift-numerics", 
    branch: "Quaternions"))

var targetDependencies: [Target.Dependency] = []
targetDependencies.append(
  .product(
    name: "Numerics",
    package: "swift-numerics"))

// MARK: - Package Manifest

targets.append(
  .target(
    name: "Mechanosynthesis",
    dependencies: targetDependencies,
    cSettings: cSettings,
    swiftSettings: swiftSettings))
targets.append(
  .testTarget(
    name: "MechanosynthesisTests",
    dependencies: ["Mechanosynthesis"],
    cSettings: cSettings,
    swiftSettings: swiftSettings))

let package = Package(
  name: "Mechanosynthesis",
  products: [
    .library(
      name: "Mechanosynthesis",
      type: .dynamic,
      targets: ["Mechanosynthesis"]),
  ],
  dependencies: packageDependencies,
  targets: targets
)
