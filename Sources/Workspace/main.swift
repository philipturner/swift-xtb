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

// Repeat this process for three slightly different singlepoints, verify that
// the results change accordingly. In addition, that the results are
// reproduced when you enter the same configuration again.
let sampleDistancesInBohr: [Float] = [
  1.6, 2.1, 3.4, 2.1
]
for distanceInBohr in sampleDistancesInBohr {
  let distanceInNm = distanceInBohr * Float(xTB_NmPerBohr)
  calculator.molecule.positions = [
    SIMD3(0.000, 0.000, 0.000),
    SIMD3(distanceInNm, 0.000, 0.000),
  ]
  
  print()
  print("basis size:", calculator.orbitals.count)
  print("occupations at finite electronic temperature:", calculator.orbitals.occupations)
  print("eigenvalues:", calculator.orbitals.eigenvalues)
  do {
    let eigenvalues = calculator.orbitals.eigenvalues
    let occupations = calculator.orbitals.occupations
    let bandEnergy = zip(eigenvalues, occupations).map(*).reduce(0, +)
    print("band energy:", bandEnergy, "zJ")
  }
  print(
    "total molecular energy:", Float(calculator.energy), "zJ", terminator: " ")
  print(
    Float(calculator.energy * xTB_HartreePerZJ), "Ha")
  
  print()
  print("distance:", distanceInNm, "nm")
  print("charges:", calculator.molecule.charges)
  print("wiberg bond orders:", calculator.molecule.bondOrders)
  print("orbital coefficents:")
  for rowID in 0..<calculator.orbitals.count {
    print("[", terminator: "")
    for columnID in 0..<calculator.orbitals.count {
      let address = rowID * calculator.orbitals.count + columnID
      let coefficient = calculator.orbitals.coefficients[address]
      var repr = String(format: "%.3f", coefficient)
      if !repr.starts(with: "-") {
        repr = " " + repr
      }
      
      if columnID == calculator.orbitals.count - 1 {
        print(repr, terminator: "")
      } else {
        print(repr, terminator: ", ")
      }
    }
    print("]")
  }
  print("forces:")
  for force in calculator.molecule.forces {
    print(force.x, "pN", terminator: " ")
    
    let scaleFactor = Float(xTB_HartreePerZJ / xTB_BohrPerNm)
    print(force.x * scaleFactor, "Ha/Bohr")
  }
}
