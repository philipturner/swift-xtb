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
  
  // The minimum number of fragments to split each electron into.
  var fragmentCount: Int?
  
  // The nuclear position for the initial guess.
  var nucleusPosition: SIMD3<Float>?
  
  // An octree initialized with the world bounds.
  var octree: Octree?
}

public struct WaveFunction {
  /// The values of the wavefunction in each cell. The cell is subdivided into
  /// a 2x2x2 group of sub-cells.
  public var cellValues: [SIMD8<Float>] = []
  
  /// The number of fragments the wavefunction should attempt to remain at.
  public var fragmentCount: Int
  
  /// The octree that stores the structure of the wavefunction.
  public var octree: Octree
  
  init(descriptor: WaveFunctionDescriptor) {
    self.fragmentCount = descriptor.fragmentCount!
    self.octree = descriptor.octree!
    
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
  }
}
