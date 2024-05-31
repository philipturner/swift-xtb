//
//  main.swift
//
//
//  Created by Philip Turner on 5/25/24.
//

import Foundation
import xTB

// Load the xTB library with a configuration for maximum performance.
setenv("OMP_STACKSIZE", "2G", 1)
setenv("OMP_NUM_THREADS", "8", 1)
xTB_Library.useLibrary(
  at: "/Users/philipturner/Documents/bypass_dependencies/libxtb.6.dylib")
try! xTB_Library.loadLibrary()

// Mute the output to the console.
xTB_Environment.verbosity = .muted

// Create a calculator.
var calculatorDesc = xTB_CalculatorDescriptor()
calculatorDesc.atomicNumbers = [7, 7]
calculatorDesc.hamiltonian = .tightBinding
let calculator = xTB_Calculator(descriptor: calculatorDesc)
calculator.molecule.positions = [
  SIMD3(0.000, 0.000, 0.000),
  SIMD3(0.110, 0.000, 0.000),
]

// Executes the singlepoint lazily, whenever a property 
// depending on it is requested.
print("basis size:", calculator.orbitals.count)
print("band energy:", 2 * calculator.orbitals.eigenvalues.reduce(0, +), "zJ")
print("total molecular energy:", calculator.energy, "zJ")
print("wiberg bond orders:", calculator.molecule.bondOrders)
print("occupations at finite electronic temperature:", calculator.orbitals.occupations)
