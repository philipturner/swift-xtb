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
// System | Homebrew OpenBLAS | Homebrew Accelerate | Custom Build |
// ------ | ----------------- | ------------------- | ------------ |
// 122    |
// 222    |
// 233    |

// Energy - GFN2-xTB
//
// System | Homebrew OpenBLAS | Homebrew Accelerate | Custom Build |
// ------ | ----------------- | ------------------- | ------------ |
// 122    |
// 222    |
// 233    |

// Latency - GFN-FF
//
// System | Homebrew OpenBLAS | Homebrew Accelerate | Custom Build |
// ------ | ----------------- | ------------------- | ------------ |
// 122    |
// 222    |
// 233    |

// Latency - GFN2-xTB
//
// System | Homebrew OpenBLAS | Homebrew Accelerate | Custom Build |
// ------ | ----------------- | ------------------- | ------------ |
// 122    |
// 222    |
// 233    |

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
xTB_Environment.verbosity = .full

// Select the system.
let system: [SIMD4<Float>] = diamondSystem233

// Create the calculator.
var calculatorDesc = xTB_CalculatorDescriptor()
calculatorDesc.atomicNumbers = system.map { UInt8($0.w) }
calculatorDesc.positions = system.map {
  SIMD3($0.x, $0.y, $0.z)
}
calculatorDesc.hamiltonian = .forceField
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
  print()
  print("actual latency:", latency)
  print("energy:", energy, "zJ")
  print("orbitals:", calculator.orbitals.count, calculator.orbitals.eigenvalues.count)
}
