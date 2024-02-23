//
//  WaveFunction.swift
//
//
//  Created by Philip Turner on 2/13/24.
//

struct WaveFunctionDescriptor {
  // The functional form of the initial guess.
  var atomicOrbital: AtomicOrbital?
  
  // The minimum number of fragments to split each electron into.
  var fragmentCount: Int?
  
  // The nuclear position for the initial guess.
  var nucleusPosition: SIMD3<Float>?
  
  // The power-2 size of the coarsest octree level.
  var sizeExponent: Int?
}

public struct WaveFunction {
  /// The values of the wavefunction at each octree node.
  public var cellValues: [SIMD8<Float>] = []
  
  /// The number of fragments the wavefunction should attempt to remain at.
  public var fragmentCount: Int
  
  /// The octree that stores the structure of the wavefunction.
  public var octree: Octree
  
  init(descriptor: WaveFunctionDescriptor) {
    guard let atomicOrbital = descriptor.atomicOrbital,
          let fragmentCount = descriptor.fragmentCount,
          let nucleusPosition = descriptor.nucleusPosition,
          let sizeExponent = descriptor.sizeExponent else {
      fatalError("Descriptor was invalid.")
    }
    self.fragmentCount = descriptor.fragmentCount!
    
    var octreeDesc = OctreeDescriptor()
    octreeDesc.sizeExponent = descriptor.sizeExponent
    self.octree = Octree(descriptor: octreeDesc)
    
    // Cache the integrals so they don't need to be recomputed. Defer the
    // initializaton of 'cellValues' until after the octree is done. That
    // roughly halves the overhead of incremental updates. Later on, we can
    // profile whether caching the wavefunctions improves performance.
    var cellIntegrals: [SIMD8<Float>] = []
    
    func createCellValues(
      center: SIMD3<Float>, spacing: Float
    ) -> SIMD8<Float> {
      var x = SIMD8<Float>(0, 1, 0, 1, 0, 1, 0, 1) * 0.5 - 0.25
      var y = SIMD8<Float>(0, 0, 1, 1, 0, 0, 1, 1) * 0.5 - 0.25
      var z = SIMD8<Float>(0, 0, 0, 0, 1, 1, 1, 1) * 0.5 - 0.25
      let centerDelta = center - nucleusPosition
      x = x * spacing + centerDelta.x
      y = y * spacing + centerDelta.y
      z = z * spacing + centerDelta.z
      
      let output = atomicOrbital.waveFunction(x: x, y: y, z: z)
      return output
    }
    
    // Perform one iteration of octree resizing.
    // - returns: Whether the octree has converged.
    func resizeOctreeNodes(_ probabilityMultiplier: Float) -> Bool {
      return false
    }
    
    // If we split at the maximum probability, the number of cells averages out
    // to roughly 4x the intended number. Therefore, we start at a looser
    // tolerance and repeat if the fragment count undershoots.
    let probabilityMultipliers: [Float] = [
      4, 2.828, 2, 1.414, 1
    ]
    for probabilityMultiplier in probabilityMultipliers {
      var iterationID = 0
      while resizeOctreeNodes(probabilityMultiplier) == false {
        iterationID += 1
        if iterationID > 50 {
          fatalError("Wave function failed to converge after 50 iterations.")
        }
      }
      
      var octreeFragmentCount = 0
      for node in octree.nodes {
        let leafMask = ~node.branchesMask
        octreeFragmentCount += leafMask.nonzeroBitCount
      }
      if octreeFragmentCount >= fragmentCount {
        
      } else if probabilityMultiplier == 1 {
        fatalError("Could not create octree with the specified fragment count.")
      }
    }
    
    #if false
    func createCellValues(metadata: SIMD4<Float>) -> SIMD8<Float> {
      var x = SIMD8<Float>(0, 1, 0, 1, 0, 1, 0, 1) * 0.5 - 0.25
      var y = SIMD8<Float>(0, 0, 1, 1, 0, 0, 1, 1) * 0.5 - 0.25
      var z = SIMD8<Float>(0, 0, 0, 0, 1, 1, 1, 1) * 0.5 - 0.25
      x = x * metadata.w + metadata.x - descriptor.nucleusPosition!.x
      y = y * metadata.w + metadata.y - descriptor.nucleusPosition!.y
      z = z * metadata.w + metadata.z - descriptor.nucleusPosition!.z
      
      let output = descriptor.atomicOrbital!.waveFunction(
        x: x, y: y, z: z)
      return output
    }
    
    // Adds the cell values for everything in the octree. Except the cells that
    // have children. Set their values to NAN.
    func fillOctree(mappedPositions: [UInt32]) {
      let NAN = SIMD8<Float>(repeating: .nan)
      var newCellValues = Array(repeating: NAN, count: octree.linkedList.count)
      for (originalPosition, mappedPosition) in mappedPositions.enumerated() {
        if mappedPosition != .max {
          newCellValues[Int(mappedPosition)] = cellValues[originalPosition]
        }
      }
      
      for nodeID in octree.linkedList.indices {
        let element = octree.linkedList[nodeID]
        if element.childCount == 0,
           newCellValues[nodeID][0].isNaN {
          let metadata = octree.metadata[nodeID]
          newCellValues[nodeID] = createCellValues(metadata: metadata)
        }
      }
      cellValues = newCellValues
    }
    
    func findCellsThatNeedChanging(probabilityMultiplier: Float) -> (
      expand: [UInt32], contract: [UInt32]
    ) {
      func createGradient(
        metadata: SIMD4<Float>, values: SIMD8<Float>
      ) -> SIMD3<Float> {
        let lowX = values.evenHalf
        let highX = values.oddHalf
        let lowZ = values.lowHalf
        let highZ = values.highHalf
        
        let casted = unsafeBitCast(values, to: SIMD4<UInt64>.self)
        let lowY = unsafeBitCast(casted.evenHalf, to: SIMD4<Float>.self)
        let highY = unsafeBitCast(casted.oddHalf, to: SIMD4<Float>.self)
        
        let deltaX = (highX - lowX).sum() / 4
        let deltaY = (highY - lowY).sum() / 4
        let deltaZ = (highZ - lowZ).sum() / 4
        let h = metadata.w
        return SIMD3(deltaX, deltaY, deltaZ) / h
      }
      
      var integralDensity: Float = .zero
      var integralGradientNorm: Float = .zero
      var cellsThatNeedExpansion: [UInt32] = []
      var cellsThatNeedContraction: Set<UInt32> = []
      var cellsThatCanContract: [UInt32] = []
      
      // First pass: calculate normalization factors.
      // Second pass: operate on normalized importance metrics.
      let maximumProbability = 1 / Float(descriptor.fragmentCount!)
      for passID in 0..<2 {
        for nodeID in octree.linkedList.indices {
          if octree.linkedList[nodeID].childCount > 0 {
            continue
          }
          
          let metadata = octree.metadata[nodeID]
          let values = cellValues[nodeID]
          let density = (values * values).sum() / 8
          let gradient = createGradient(metadata: metadata, values: values)
          let gradientNorm = (gradient * gradient).sum()
          let volume = metadata.w * metadata.w * metadata.w
          
          let ε: Float = .leastNormalMagnitude
          let densityContribution = max(ε, density * volume)
          let gradientContribution = max(ε, gradientNorm * volume)
          if passID == 0 {
            integralDensity += densityContribution
            integralGradientNorm += gradientContribution
          } else {
            var importanceMetric: Float = .zero
            importanceMetric += densityContribution / integralDensity
            importanceMetric += gradientContribution / integralGradientNorm
            importanceMetric *= 0.5
            
            let upperBound = probabilityMultiplier * maximumProbability
            let lowerBound = probabilityMultiplier * maximumProbability / 8
            if importanceMetric > upperBound {
              cellsThatNeedExpansion.append(UInt32(nodeID))
            } else if importanceMetric < lowerBound {
              cellsThatNeedContraction.insert(UInt32(nodeID))
            }
          }
        }
      }
      
    outer:
      for nodeID in octree.linkedList.indices {
        if octree.linkedList[nodeID].childCount == 0 {
          continue
        }
        
        for i in 1...8 {
          if octree.linkedList[nodeID + i].childCount == 8 {
            continue outer
          }
          let childID = UInt32(nodeID + i)
          guard cellsThatNeedContraction.contains(childID) else {
            continue outer
          }
        }
        cellsThatCanContract.append(UInt32(nodeID))
      }
      
      return (cellsThatNeedExpansion, cellsThatCanContract)
    }
    
    fillOctree(mappedPositions: [])
    
    // If we split at the maximum probability, the number of cells averages out
    // to roughly 4x the intended number. Therefore, we start at a looser
    // tolerance and repeat if the fragment count undershoots.
    let probabilityMultipliers: [Float] = [
      4, 2.828, 2, 1.414, 1
    ]
    for probabilityMultiplier in probabilityMultipliers {
      var converged = false
      for _ in 0..<100 {
        let (expand, contract) = findCellsThatNeedChanging(
          probabilityMultiplier: probabilityMultiplier)
        if expand.isEmpty, contract.isEmpty {
          converged = true
          break
        }
        
        let mappedPositions = octree.modifyNodes(
          expand: expand, contract: contract)
        fillOctree(mappedPositions: mappedPositions)
      }
      guard converged else {
        fatalError("Wave function failed to converge after 100 iterations.")
      }
      
      // WARNING: This is not correct! The parent nodes should be excluded from
      // the count. In addition, the 2x2x2 sub-cells should be reported to most
      // accurately reflect the compute cost. Revise the test suite to account
      // for the new definition of 'fragment'. This includes increasing the
      // high-accuracy measurements from 1,000,000 -> 8,000,000 fragments.
      if octree.linkedList.count >= descriptor.fragmentCount! {
        break
      }
    }
    #endif
    
    fatalError("Not implemented.")
  }
}
