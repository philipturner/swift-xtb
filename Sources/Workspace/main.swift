//
//  main.swift
//
//
//  Created by Philip Turner on 5/25/24.
//

import Darwin
import xTB

// Load the xTB library with a configuration for maximum performance.
setenv("OMP_STACKSIZE", "2G", 1)
setenv("OMP_NUM_THREADS", "8", 1)
xTB_Library.useLibrary(
  at: "/Users/philipturner/Documents/OpenMM/bypass_dependencies/libxtb.6.dylib")
try! xTB_Library.loadLibrary()

// Create an environment.
xTB_Environment.verbosity = .muted
// xTB_Environment.show()

// Create a calculator.
var calculatorDesc = xTB_CalculatorDescriptor()
calculatorDesc.atomicNumbers = [7, 7]
// calculatorDesc.positions = [
//   SIMD3(0.000, 0.000, 0.000),
//   SIMD3(0.110, 0.000, 0.000),
// ]
let calculator = xTB_Calculator(descriptor: calculatorDesc)

// Set the external charges.
calculator.externalCharges.atomicNumbers = []
calculator.externalCharges.charges = []
calculator.externalCharges.positions = []
for chargeID in 0..<2 {
  let atomicNumber: UInt8 = 7
  let charge: Float = (chargeID == 0) ? 1 : -1
  var position: SIMD3<Float>
  if chargeID == 0 {
    position = SIMD3(0.220, 0.000, 0.000)
  } else {
    position = SIMD3(-0.110, 0.000, 0.000)
  }
  
  calculator.externalCharges.atomicNumbers.append(atomicNumber)
  calculator.externalCharges.charges.append(charge)
  calculator.externalCharges.positions.append(position)
}
calculator.molecule.positions = [
  SIMD3(0.000, 0.000, 0.000),
  SIMD3(0.111, 0.000, 0.000),
]

// Run a singlepoint.
calculator.temporaryTestFunction()
xTB_Environment.show()
