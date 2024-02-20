//
//  WaveFunction.swift
//
//
//  Created by Philip Turner on 2/13/24.
//

// Stores a set of fragments. Each fragment internally stores a 2x2x2 group
// of sub-fragments contiguously. One can compute the first derivative at the
// center of the cell for XC functionals.

struct WaveFunctionDescriptor {
  // The functional form of the initial guess.
  var atomicOrbital: AtomicOrbital?
  
  // The nuclear position for the initial guess.
  var nucleusPosition: SIMD3<Float>?
  
  // The minimum number of fragments to split the wavefunction into.
  //
  // The heuristic for determining the importance of a specific 3D location is
  // not specified. It may depend on charge density or wavefunction gradient.
  var minimumFragmentCount: Float?
  
  // The data to initialize the octree descriptor with.
  var worldBounds: SIMD4<Float>?
}

struct WaveFunction {
  // The values of the wavefunction in each cell. The cell is subdivided into
  // a 2x2x2 group of sub-cells.
  var cellValues: [SIMD8<Float>] = []
  
  // The octree that stores the structure of the wavefunction.
  var octree: Octree
  
  init(descriptor: WaveFunctionDescriptor) {
    var octreeDesc = OctreeDescriptor()
    octreeDesc.origin = unsafeBitCast(descriptor.worldBounds!, to: SIMD3.self)
    octreeDesc.size = descriptor.worldBounds!.w
    self.octree = Octree(descriptor: octreeDesc)
    
    // TODO: Start determining how 'Octree' will be used, by filling one in
    // with the initial guess to an orbital. The very first node is centered
    // at the origin. Create a recursive procedure that keeps adding nodes. It
    // terminates once the constraints for importance metrics are satisfied.
    
    func createCellValues(metadata: SIMD4<Float>) -> SIMD8<Float> {
      let origin = unsafeBitCast(metadata, to: SIMD3<Float>.self)
      let size = metadata.w
      var output: SIMD8<Float> = .zero
      
      for childID in 0..<8 {
        let xIndex = UInt32(childID) % 2
        let yIndex = UInt32(childID >> 1) % 2
        let zIndex = UInt32(childID >> 2) % 2
        var delta: SIMD3<Float> = .init(repeating: -0.25)
        let indices = SIMD3<UInt32>(xIndex, yIndex, zIndex)
        delta.replace(with: 0.25, where: indices .> 0)
        
        let position = origin + delta * size
        let nucleusDelta = position - descriptor.nucleusPosition!
        let waveFunction = descriptor.atomicOrbital!
          .waveFunction(position: nucleusDelta)
        output[childID] = waveFunction
      }
      return output
    }
    
    // Adds the cell values for everything in the octree. Except the cells that
    // have children. Their values are set to NAN.
    func fillOctree() {
      cellValues = []
      for nodeID in octree.linkedList.indices {
        let element = octree.linkedList[nodeID]
        var newValue: SIMD8<Float>
        
        if element.childCount == 8 {
          newValue = .init(repeating: .nan)
        } else {
          let metadata = octree.metadata[nodeID]
          newValue = createCellValues(metadata: metadata)
        }
        cellValues.append(newValue)
      }
    }
    
    func markCellsThatNeedExpansion() {
      /*
      for nodeID in octree.linkedList.indices {
        let cellValue = octree.cellValues[nodeID]
        let density = (cellValue * cellValue).sum() / 8
        let radius = Float.exp2(Float(currentHierarchyLevel))
        let probability = density * radius * radius * radius
        let maximumProbability = 1 / descriptor.minimumFragmentCount!
        
        if probability > maximumProbability {
          octree.cellValues[nodeID] = .init(repeating: .nan)
          octree.appendChildren(to: [UInt32(nodeID)])
          // The children don't have cell values now.
        }
      }
       */
    }
    
    func markCellsThatNeedContraction() {
      
    }
  }
}
