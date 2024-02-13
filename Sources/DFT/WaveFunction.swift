//
//  WaveFunction.swift
//
//
//  Created by Philip Turner on 2/13/24.
//

import Foundation

// Stores a set of fragments. Each fragment internally stores a 2x2x2 group
// of sub-fragments contiguously. It can compute the first derivative at the
// center of the cell for XC functionals. It can also compute the second
// derivative at each sub-fragment, provided information about nearby fragments.

struct WaveFunctionDescriptor {
  // The functional form of the initial guess.
  var atomicOrbital: AtomicOrbital?
  
  // The nuclear position for the initial guess.
  var position: SIMD3<Float>?
  
  // The maximum probability contained by any fragment.
  var maximumProbability: Float?
  
  // The diameter of the highest multigrid level. For now, bounding boxes must
  // be power-2 cubes centered at the origin.
  var diameter: Float?
}

struct WaveFunction {
  // Array of wave function values. The multigrid stores indices that point to
  // the array index of the 2x2x2 supercell.
  var values: [SIMD8<Float>] = []
  
  init(descriptor: WaveFunctionDescriptor) {
    guard let atomicOrbital = descriptor.atomicOrbital,
          let position = descriptor.position,
          let maximumProbability = descriptor.maximumProbability,
          let diameter = descriptor.diameter else {
      fatalError("Wave function descriptor not complete.")
    }
  }
}
