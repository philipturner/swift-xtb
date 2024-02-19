//
//  WaveFunction.swift
//
//
//  Created by Philip Turner on 2/13/24.
//

struct WaveFunctionDescriptor {
  // The functional form of the initial guess.
  var atomicOrbital: AtomicOrbital?
  
  // The radius of the highest octree level, in powers of 2. The octree is
  // centered at the origin.
  var highestLevelSize: Int?
  
  // The nuclear position for the initial guess.
  var position: SIMD3<Float>?
  
  // The minimum number of fragments to split the wavefunction into.
  //
  // The heuristic for determining the importance of a specific 3D location is
  // not specified. It may depend on charge density or wavefunction gradient.
  var minimumFragmentCount: Float?
}

struct WaveFunction {
  var octree: Octree
  
  init(descriptor: WaveFunctionDescriptor) {
    self.octree = Octree()
    
    // TODO: Start determining how 'Octree' will be used, by filling one in
    // with the initial guess to an orbital. The very first node is centered
    // at the origin. Create a recursive procedure that keeps adding nodes. It
    // terminates once the constraints for importance metrics are satisfied.
    //
    // For the first sketch, make the importance metric be un-normalized
    // density for simplicity. After you've got all the code figured out, go
    // back and replace with normalized density + gradient^2.
    
    for childID in 0..<8 {
      let xIndex = UInt32(childID) % 2
      let yIndex = UInt32(childID >> 1) % 2
      let zIndex = UInt32(childID >> 2) % 2
      let radius = Float.exp2(Float(descriptor.highestLevelSize!))
      let lowerCorner = SIMD3<Float>(repeating: -radius)
      
      var position = lowerCorner / 2
      var delta: SIMD3<Float> = .zero
      let indices = SIMD3<UInt32>(xIndex, yIndex, zIndex)
      delta.replace(with: .one, where: indices .> 0)
      position += delta * radius
      
      let nucleusDelta = position - descriptor.position!
      let waveFunction = descriptor.atomicOrbital!
        .waveFunction(position: nucleusDelta)
      octree.cellValues[0][childID] = waveFunction
    }
    
    
  }
}
