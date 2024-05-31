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

// Create a calculator.
var calculatorDesc = xTB_CalculatorDescriptor()
calculatorDesc.atomicNumbers = [7, 7]
calculatorDesc.positions = [
  SIMD3(0.000, 0.000, 0.000),
  SIMD3(0.111, 0.000, 0.000),
]
calculatorDesc.hamiltonian = .tightBinding
let calculator = xTB_Calculator(descriptor: calculatorDesc)

// TODO: Get the Wiberg bond orders for the N2 molecule.
//
// TODO: Plot out the potential energy surface for N2 across various bond
// lengths, compare to the graph from INQ.
xTB_Environment.show()
