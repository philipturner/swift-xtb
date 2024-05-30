//
//  main.swift
//
//
//  Created by Philip Turner on 5/25/24.
//

import Darwin
import xTB

// From the command line:
//   swift run -Xswiftc -Ounchecked --debugger Workspace -o run
// Alternatively, change the Xcode target to the 'Workspace' executable.

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
