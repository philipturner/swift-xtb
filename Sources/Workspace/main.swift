//
//  main.swift
//
//
//  Created by Philip Turner on 5/25/24.
//

import Foundation
import Numerics
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
  SIMD3(0.180, 0.000, 0.000),
]

// Make a constant for atom count, for convenience.
let atomCount = calculator.molecule.atomicNumbers.count

// FIRE Algorithm
var Δt: Float = 0.001
var NP0: Int = 0
var oldState: [SIMD3<Float>]?
var velocities = [SIMD3<Float>](repeating: .zero, count: atomCount)

// Loop until the maximum number of iterations is reached.
for frameID in 0..<10 {
  defer { print() }
  print("frame: \(frameID)", terminator: " | ")
  print("energy:", Float(calculator.energy), "zJ", terminator: " | ")
  print("Δt:", Δt, terminator: " | ")
  
  // Find the power (P) and maximum force.
  let forces = calculator.molecule.forces
  var P: Float = .zero
  var maxForce: Float = .zero
  for atomID in 0..<atomCount {
    let force = forces[atomID]
    let velocity = velocities[atomID]
    P += (force * velocity).sum()
    
    let forceMagnitude = (force * force).sum().squareRoot()
    maxForce = max(maxForce, forceMagnitude)
  }
  print("max force:", maxForce, terminator: " | ")
  
  // Either restart or increase the timestep.
  if frameID > 0, P < 0 {
    if let oldState {
      calculator.molecule.positions = oldState
    }
    velocities = [SIMD3<Float>](repeating: .zero, count: atomCount)
    
    NP0 = 0
    Δt = max(0.5 * Δt, 0.001 / 20)
  } else {
    NP0 += 1
    if NP0 > 5 {
      Δt = min(1.1 * Δt, 0.005)
    }
  }
  
  // Update the atoms.
  oldState = []
  for atomID in 0..<atomCount {
    let force = forces[atomID]
    var velocity = velocities[atomID]
    
    // Ensure the force scale is finite.
    var forceScale: Float
    do {
      let vNorm = (velocity * velocity).sum().squareRoot()
      let fNorm = (force * force).sum().squareRoot()
      forceScale = vNorm / fNorm
    }
    if forceScale.isNaN || forceScale.isInfinite {
      forceScale = .zero
    }
    
    // Choose a mass for the atom.
    var mass: Float
    do {
      let atomicNumber = calculator.molecule.atomicNumbers[atomID]
      if atomicNumber == 1 {
        // Set the hydrogen mass to 4 amu.
        mass = 4 * 1.6605
      } else {
        // Normalize the vibration period, to match carbon.
        mass = 12.011 * 1.6605
      }
    }
    
    // Semi-implicit Euler integration.
    let α: Float = 0.1
    velocity += Δt * force / mass
    velocity = (1 - α) * velocity + α * force * forceScale
    
    // Integrate the position.
    let position = calculator.molecule.positions[atomID]
    let halfwayPoint = position + 0.5 * Δt * velocity
    let newPosition = position + Δt * velocity
    print(velocity.x, velocity.y, velocity.z, terminator: " | ")
    
    // Store the new state.
    velocities[atomID] = velocity
    calculator.molecule.positions[atomID] = newPosition
    oldState!.append(halfwayPoint)
    
    let delta = Δt * velocity
    let distance = (delta * delta).sum().squareRoot()
    print("d =", distance, terminator: " | ")
  }
}

// Report the final state.
let positions = calculator.molecule.positions
print("final bond vector:", positions[1] - positions[0])
