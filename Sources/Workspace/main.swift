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

// Try out the xTB_Environment API.
let environment = xTB_Environment()
environment.verbosity = .full
environment.show()

// Try out the xTB_Calculator API.
var calculatorDesc = xTB_CalculatorDescriptor()
calculatorDesc.atomicNumbers = [7, 7]
calculatorDesc.environment = environment

// Create a calculator.
let calculator = xTB_Calculator(descriptor: calculatorDesc)
calculator.setPositions([
  SIMD3(0.000, 0.000, 0.000),
  SIMD3(0.110, 0.000, 0.000),
])

// Try out the xTB_ExternalCharge API.
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
}
calculator.setExternalCharges(externalCharges)
environment.show()

/*

// Try out the xTB_Molecule API.
var moleculeDesc = xTB_MoleculeDescriptor()
moleculeDesc.atomicNumbers = [7, 7]
moleculeDesc.environment = environment

// Create a molecule.
let molecule = xTB_Molecule(descriptor: moleculeDesc)
molecule.setPositions([
  SIMD3(0.000, 0.000, 0.000),
  SIMD3(0.110, 0.000, 0.000),
])
print(environment.status)

// Create a calculator.
var calculatorDesc = xTB_CalculatorDescriptor()
calculatorDesc.molecule = molecule
let calculator = xTB_Calculator(descriptor: calculatorDesc)
print(environment.status)

*/
