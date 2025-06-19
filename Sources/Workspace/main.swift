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

// WARNING: Watch out for 'gfnff_topo' files leaking into the working directory.
// Perhaps the file doesn't appear when you set verbosity to muted?
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
