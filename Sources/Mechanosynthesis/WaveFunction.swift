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
  
  /// The octree that stores the structure of the wavefunction.
  public var octree: Octree
  
  init(descriptor: WaveFunctionDescriptor) {
    self.octree = descriptor.octree!
    
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
    
    // There isn't a major need to contract too-small cells during
    // initialization. The gradual increase in normalization factor only
    // slightly lowers the importance value.
    func findCellsThatNeedExpansion(probabilityMultiplier: Float) -> [UInt32] {
      // This might eventually be brought out into a global function.
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
      var output: [UInt32] = []
      
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
            
            if importanceMetric > probabilityMultiplier * maximumProbability {
              output.append(UInt32(nodeID))
            }
          }
        }
      }
      
      return output
    }
    
    fillOctree()
    
    // If we split at the maximum probability, the number of cells
    // averages out to roughly 4x the intended number. Therefore, we start at
    // a looser tolerance and repeat if the fragment count undershoots.
    let probabilityMultipliers: [Float] = [
      5.656, 4, 2.828, 2, 1.414, 1
    ]
    for probabilityMultiplier in probabilityMultipliers {
      var converged = false
      for _ in 0..<100 {
        let cellsToExpand = findCellsThatNeedExpansion(
          probabilityMultiplier: probabilityMultiplier)
        if cellsToExpand.isEmpty {
          converged = true
          break
        }
        octree.insertChildren(at: cellsToExpand)
        fillOctree()
      }
      guard converged else {
        fatalError("Wave function failed to converge after 100 iterations.")
      }
      if octree.linkedList.count >= descriptor.fragmentCount! {
        break
      }
    }
  }
}
