//
//  main.swift
//  xTB
//
//  Created by Philip Turner on 6/19/25.
//

import Foundation
import xTB

// Try making a fork of homebrew-qc and installing it on my system. Updated
// with mctc v0.4.0, which may fix issues causing the crash.



// Energy - GFN-FF
//
// System | Homebrew OpenBLAS | Homebrew Accelerate | Latest Commit    |
// ------ | ----------------- | ------------------- | ---------------- |
// 122    |     -57126.394 zJ |       -57126.394 zJ |    -57126.394 zJ |
// 222    |    -100762.569 zJ |      -100762.569 zJ |   -100762.569 zJ |
// 233    |    -208433.751 zJ |      -208433.751 zJ |   -208433.751 zJ |

// Latency - GFN-FF
//
// System | Homebrew OpenBLAS | Homebrew Accelerate | Latest Commit    |
// ------ | ----------------- | ------------------- | ---------------- |
// 122    |            1.6 ms |              1.5 ms |           1.2 ms |
// 222    |            2.9 ms |              3.0 ms |           2.4 ms |
// 233    |           10.0 ms |              9.6 ms |           8.2 ms |



// Energy - GFN2-xTB
//
// System | Homebrew OpenBLAS | Homebrew Accelerate | Latest Commit    |
// ------ | ----------------- | ------------------- | ---------------- |
// 122    |    -452321.592 zJ |      -452321.592 zJ |   -452321.592 zJ |
// 222    |    -841344.658 zJ |      -841344.658 zJ |   -841344.658 zJ |
// 233    |           crashes |             crashes |          crashes |

// Latency - GFN2-xTB
//
// System | Homebrew OpenBLAS | Homebrew Accelerate | Latest Commit    |
// ------ | ----------------- | ------------------- | ---------------- |
// 122    |          268.3 ms |             70.5 ms |          63.2 ms |
// 222    |          870.1 ms |            264.8 ms |         219.4 ms |
// 233    |           crashes |             crashes |          crashes |



// Prepare the environment for maximum performance with xTB.
setenv("OMP_STACKSIZE", "2G", 1)
setenv("OMP_NUM_THREADS", "8", 1) // replace '8' with number of P-cores

// Suppress unwanted output from GFN-FF.
let url = FileManager.default.temporaryDirectory
let path = url.relativePath
let worked = FileManager.default.changeCurrentDirectoryPath(path)
guard worked else {
  fatalError("Could not redirect gfnff_topo directory.")
}
xTB_Environment.verbosity = .muted
xTB_Environment.setOutput("/dev/null")

// Select the system.
let system: [SIMD4<Float>] = diamondSystem222

// Create the calculator.
var calculatorDesc = xTB_CalculatorDescriptor()
calculatorDesc.atomicNumbers = system.map { UInt8($0.w) }
calculatorDesc.positions = system.map {
  SIMD3($0.x, $0.y, $0.z)
}
//calculatorDesc.hamiltonian = .forceField
let calculator = xTB_Calculator(descriptor: calculatorDesc)

// Run just one loop iteration.
for _ in 0..<1 {
  calculator.molecule.positions = system.map {
    SIMD3($0.x, $0.y, $0.z)
  }
  
  let checkpoint0 = Date()
  let energy = calculator.energy
  let checkpoint1 = Date()
  let latency = checkpoint1.timeIntervalSince(checkpoint0)
  
  let formattedLatency = String(format: "%.1f", latency * 1000)
  print()
  print("actual latency:", formattedLatency, "ms")
  
  let formattedEnergy = String(format: "%.3f", energy)
  print("energy:", formattedEnergy, "zJ")
  print("orbitals:", calculator.orbitals.count, calculator.orbitals.eigenvalues.count)
}
