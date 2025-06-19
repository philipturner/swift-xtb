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

/*
// Load the 'xtb' dylib.
xTB_Library.useLibrary(
  at: "/Users/philipturner/Documents/MolecularRenderer/swift-xtb/libxtb.6.dylib")
try! xTB_Library.loadLibrary()
*/

// WARNING: Watch out for 'gfnff_topo' files leaking into the working directory.
// Perhaps the file doesn't appear when you set verbosity to muted?
xTB_Environment.verbosity = .minimal

// Create the calculator.
var calculatorDesc = xTB_CalculatorDescriptor()
calculatorDesc.atomicNumbers = diamondSystem122.map { UInt8($0.w) }
calculatorDesc.positions = diamondSystem122.map {
  SIMD3($0.x, $0.y, $0.z)
}
let calculator = xTB_Calculator(descriptor: calculatorDesc)

// Run just one loop iteration.
for _ in 0..<1 {
  calculator.molecule.positions = diamondSystem122.map {
    SIMD3($0.x, $0.y, $0.z)
  }
  
  let checkpoint0 = Date()
  _ = calculator.energy
  let checkpoint1 = Date()
  let latency = checkpoint1.timeIntervalSince(checkpoint0)
  print()
  print("actual latency:", latency)
  print("energy:", calculator.energy)
}
