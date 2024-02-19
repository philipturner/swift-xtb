//
//  WaveFunction.swift
//
//
//  Created by Philip Turner on 2/13/24.
//

struct WaveFunctionDescriptor {
  // The functional form of the initial guess.
  var atomicOrbital: AtomicOrbital?
  
  
  
  // The nuclear position for the initial guess.
  var position: SIMD3<Float>?
  
  // The minimum number of fragments to split the wavefunction into.
  //
  // The heuristic for determining the importance of a specific 3D location is
  // not specified. It may depend on charge density or wavefunction gradient.
  var minimumFragmentCount: Float?
}

struct WaveFunction {
  // The radius of the highest octree level, in powers of 2. The octree is
  // centered at the origin.
  var highestLevelSize: Int
  
  init(descriptor: WaveFunctionDescriptor) {
    // Generate an octree and the appropriate bounds from the atomic orbital.
    // The bounds may change throughout the self-consistent cycle.
    fatalError("Not implemented.")
  }
}
