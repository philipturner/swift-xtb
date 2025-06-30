// swift-tools-version: 6.1

import PackageDescription
import class Foundation.ProcessInfo

var linkerSettings: [LinkerSetting] = []

if let path = ProcessInfo.processInfo.environment["XTB_LIBRARY_PATH"] {
  linkerSettings += [
    .unsafeFlags(["-L\(path)"]),
    .linkedLibrary("xtb"),
  ]
  
  #if os(Windows)
  linkerSettings += [
    .linkedLibrary("libgcc_s_seh-1"),
    .linkedLibrary("libgfortran-5"),
    .linkedLibrary("libgomp-1"),
    .linkedLibrary("libopenblas"),
    .linkedLibrary("libquadmath-0"),
    .linkedLibrary("libwinpthread-1"),
    .linkedLibrary("libmingwex"),
    .linkedLibrary("libgcc"),
  ]
  #endif
}

let package = Package(
  name: "swift-xtb",
  dependencies: [],
  targets: [
    .target(
      name: "C_xTB",
      dependencies: [],
      linkerSettings: linkerSettings),
    .target(
      name: "xTB",
      dependencies: ["C_xTB"]),
    .executableTarget(
      name: "Workspace",
      dependencies: ["xTB"]),
  ]
)
