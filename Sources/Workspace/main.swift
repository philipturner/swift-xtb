//
//  main.swift
//  xTB
//
//  Created by Philip Turner on 6/19/25.
//

import Foundation
import xTB

var limit = rlimit()
let returnValue = getrlimit(RLIMIT_STACK, &limit)
guard returnValue == 0 else {
  fatalError("Could not get rlimit. Encountered error code: \(returnValue)")
}
print("cur:", limit.rlim_cur as UInt64)
print("max:", limit.rlim_max as UInt64)

/*
 8176 KB
 cur: 8_372_224
 max: 67_092_480
 */

/*
 65520 KB
 cur: 67_092_480
 max: 67_092_480
 */

// Prepare the environment for maximum performance with xTB.
setenv("OMP_STACKSIZE", "16G", 1)
setenv("OMP_NUM_THREADS", "1", 1) // replace '8' with number of P-cores

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
calculatorDesc.hamiltonian = .tightBinding
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
  print("energy:", energy)
  print("orbitals:", calculator.orbitals.count, calculator.orbitals.eigenvalues.count)
}
