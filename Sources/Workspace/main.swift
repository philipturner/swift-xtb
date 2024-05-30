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
let environment = xTB_Environment()
environment.verbosity = .full
environment.show()

// Create a calculator.
var calculatorDesc = xTB_CalculatorDescriptor()
calculatorDesc.atomicNumbers = [7, 7]
calculatorDesc.environment = environment
let calculator = xTB_Calculator(descriptor: calculatorDesc)

// Set the positions.
calculator.positions = [
  SIMD3(0.000, 0.000, 0.000),
  SIMD3(0.110, 0.000, 0.000),
]

// Set the external charges.
var externalCharges: [xTB_ExternalCharge] = []
for chargeID in 0..<2 {
  var charge = xTB_ExternalCharge()
  charge.charge = (chargeID == 0) ? 1 : -1
  charge.chemicalHardness = 99
  
  if chargeID == 0 {
    charge.position = SIMD3(0.220, 0.000, 0.000)
  } else {
    charge.position = SIMD3(-0.110, 0.000, 0.000)
  }
  externalCharges.append(charge)
}
calculator.externalCharges = externalCharges
calculator.externalCharges = externalCharges
environment.show()
