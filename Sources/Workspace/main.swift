//
//  main.swift
//  xTB
//
//  Created by Philip Turner on 6/19/25.
//

import Foundation
import xTB

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

// Create the calculator.
var calculatorDesc = xTB_CalculatorDescriptor()
calculatorDesc.atomicNumbers = diamondSystem122.map { UInt8($0.w) }
calculatorDesc.positions = diamondSystem122.map {
  SIMD3($0.x, $0.y, $0.z)
}
calculatorDesc.hamiltonian = .forceField
let calculator = xTB_Calculator(descriptor: calculatorDesc)

// Run just one loop iteration.
for _ in 0..<10 {
  calculator.molecule.positions = diamondSystem122.map {
    SIMD3($0.x, $0.y, $0.z)
  }
  
  let checkpoint0 = Date()
  _ = calculator.energy
  let checkpoint1 = Date()
  let latency = checkpoint1.timeIntervalSince(checkpoint0)
  print()
  print("actual latency:", latency)
}
