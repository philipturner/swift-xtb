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

// TODO: Plot out the potential energy surface for N2 across various bond
// lengths, compare to the graph from INQ.
for stepID in 0..<20 {
  let positionInBohr = 1.6 + 0.1 * Float(stepID)
  let positionInNm = positionInBohr * Float(xTB_NmPerBohr)
  calculator.molecule.positions = [
    SIMD3(0.000, 0.000, 0.000),
    SIMD3(positionInNm, 0.000, 0.000),
  ]
  
  let energyInHa = 5.8 + Float(calculator.energy * xTB_HartreePerZJ)
  print(String(format: "%.3f", positionInBohr), "Bohr", terminator: " | ")
  print(String(format: "%.3f", energyInHa), "Hartree")
}
xTB_Environment.show()

// TODO: Get the Wiberg bond orders for the N2 molecule.
