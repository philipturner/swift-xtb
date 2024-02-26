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
//  /// The values of the wavefunction at each octree node.
//  public var cellValues: [SIMD8<Float>] = []
  
  // Replacing cell values with the atomic orbital for now, for debugging
  // purposes.
  var atomicOrbital: AtomicOrbital
  public func atomicOrbitalWaveFunction(
    x: SIMD8<Float>,
    y: SIMD8<Float>,
    z: SIMD8<Float>
  ) -> SIMD8<Float> {
    atomicOrbital.waveFunction(x: x, y: y, z: z)
  }
  
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
    self.atomicOrbital = atomicOrbital
    self.fragmentCount = fragmentCount
    
    var octreeDesc = OctreeDescriptor()
    octreeDesc.sizeExponent = sizeExponent
    self.octree = Octree(descriptor: octreeDesc)
    
    // Cache the integrals so they don't need to be recomputed. Defer the
    // initialization of 'cellValues' until after the octree is done. We can
    // profile whether caching the wavefunctions improves performance.
    var densityIntegrals: [SIMD8<Float>] = []
    var gradientIntegrals: [SIMD8<Float>] = []
    densityIntegrals.append(SIMD8(repeating: .nan))
    gradientIntegrals.append(SIMD8(repeating: .nan))
    
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
    
    // Evaluate the density and square gradient over all of real space. Operate
    // on the caches, 'densityIntegrals' and 'gradientIntegrals'.
    func evaluateGlobalIntegrals() -> (
      density: Float,
      squareGradient: Float
    ) {
      var density: Double = .zero
      var squareGradient: Double = .zero
      for nodeID in octree.nodes.indices {
        var densityIntegral = densityIntegrals[nodeID]
        var gradientIntegral = gradientIntegrals[nodeID]
        let node = octree.nodes[nodeID]
        
        // Compute the integrals and cache them.
        if densityIntegrals[nodeID][0].isNaN {
          let Ψ = createCellValues(center: node.center, spacing: node.spacing)
          let d3r = node.spacing * node.spacing * node.spacing / 8
          densityIntegral = Ψ * 1 * Ψ * d3r
          
          let casted = unsafeBitCast(Ψ, to: SIMD4<UInt64>.self)
          let lowX = Ψ.evenHalf
          let highX = Ψ.oddHalf
          let lowY = unsafeBitCast(casted.evenHalf, to: SIMD4<Float>.self)
          let highY = unsafeBitCast(casted.oddHalf, to: SIMD4<Float>.self)
          let lowZ = Ψ.lowHalf
          let highZ = Ψ.highHalf
          
          let Δx = (highX - lowX) / node.spacing
          let Δy = (highY - lowY) / node.spacing
          let Δz = (highZ - lowZ) / node.spacing
          let intΔx = Δx * Δx * d3r
          let intΔy = Δy * Δy * d3r
          let intΔz = Δz * Δz * d3r
          
          let castedIntΔx = unsafeBitCast(intΔx, to: SIMD4<UInt32>.self)
          var castedIntΔx2 = SIMD4<UInt64>(truncatingIfNeeded: castedIntΔx)
          castedIntΔx2 |= castedIntΔx2 &<< 32
          let castedIntΔy = unsafeBitCast(intΔy, to: SIMD2<UInt64>.self)
          let castedIntΔy2 = SIMD4(castedIntΔy[0], castedIntΔy[0],
                                   castedIntΔy[1], castedIntΔy[1])
          
          let finalIntΔx = unsafeBitCast(castedIntΔx2, to: SIMD8<Float>.self)
          let finalIntΔy = unsafeBitCast(castedIntΔy2, to: SIMD8<Float>.self)
          let finalIntΔz = SIMD8<Float>(lowHalf: intΔz, highHalf: intΔz)
          gradientIntegral = finalIntΔx + finalIntΔy + finalIntΔz
          
          // Check the correctness of this optimized code by calculating it the
          // slow way.
          for laneID in 0..<8 {
            // The sign of the gradient doesn't matter; it will be squared.
            let partX = Ψ[laneID ^ 1] - Ψ[laneID]
            let partY = Ψ[laneID ^ 2] - Ψ[laneID]
            let partZ = Ψ[laneID ^ 4] - Ψ[laneID]
            var partXYZ = SIMD3(partX, partY, partZ)
            partXYZ /= node.spacing
            partXYZ = partXYZ * partXYZ * d3r
            
            let expected = partXYZ.sum()
            let actual = gradientIntegral[laneID]
            guard expected == actual else {
              let absoluteDifference = actual - expected
              let relativeDifference = actual / expected
              fatalError("Incorrect gradient integral. Expected \(expected), got \(actual) (difference: \(absoluteDifference), ratio: \(relativeDifference)).")
            }
          }
          
          let ε: Float = .leastNormalMagnitude
          densityIntegral.replace(with: ε, where: densityIntegral .< ε)
          gradientIntegral.replace(with: ε, where: gradientIntegral .< ε)
          densityIntegrals[nodeID] = densityIntegral
          gradientIntegrals[nodeID] = gradientIntegral
        }
        
        let shifts = SIMD8<UInt8>(0, 1, 2, 3, 4, 5, 6, 7)
        var mask8 = SIMD8<UInt8>(repeating: 1) &<< shifts
        mask8 = mask8 & node.branchesMask
        let mask32 = SIMD8<UInt32>(truncatingIfNeeded: mask8) .!= 0
        densityIntegral.replace(with: 0, where: mask32)
        gradientIntegral.replace(with: 0, where: mask32)
        density += Double(densityIntegral.sum())
        squareGradient += Double(gradientIntegral.sum())
      }
      return (Float(density), Float(squareGradient))
    }
    
    var consoleMessage: String = ""
    
    // Perform one iteration of octree resizing.
    func queryOctreeNodes(
      probabilityMultiplier: Float
    ) -> (
      expanded: [(UInt32, UInt8)],
      contracted: [UInt32]
    ) {
      let (globalDensity, globalSquareGradient) = evaluateGlobalIntegrals()
      var expanded: [(UInt32, UInt8)] = []
      var contracted: [UInt32] = []
      
      for nodeID in octree.nodes.indices {
        var importanceMetric: SIMD8<Float> = .zero
        importanceMetric += densityIntegrals[nodeID] / globalDensity
        importanceMetric += gradientIntegrals[nodeID] / globalSquareGradient
        importanceMetric *= 0.5
        
        let branchesMask = octree.nodes[nodeID].branchesMask
        let shifts = SIMD8<UInt8>(0, 1, 2, 3, 4, 5, 6, 7)
        let mask8 = (SIMD8<UInt8>(repeating: 1) &<< shifts) & branchesMask
        let mask32 = SIMD8<UInt32>(truncatingIfNeeded: mask8) .!= 0
        importanceMetric.replace(with: 0, where: mask32)
        
        let threshold = probabilityMultiplier / Float(fragmentCount)
        
        var thresholdMask32 = SIMD8<UInt32>(repeating: .zero)
        thresholdMask32.replace(with: 1, where: importanceMetric .> threshold)
        var thresholdMask8 = SIMD8<UInt8>(truncatingIfNeeded: thresholdMask32)
        thresholdMask8 = thresholdMask8 &<< shifts
        let thresholdMaskSum = thresholdMask8.wrappedSum()
        
        /*
         expanded 0 contracted 4232
         expanded 7 contracted 505
         expanded 1 contracted 78
         expanded 7 contracted 0
         expanded 0 contracted 21
         expanded 7 contracted 0
         expanded 0 contracted 21
         expanded 7 contracted 0
         expanded 0 contracted 21
         expanded 7 contracted 0
         probability: 4.0, fragment count: 37724, expected: 40000
         expanded 743 contracted 0
         expanded 0 contracted 0
         expanded 0 contracted 53
         expanded 39 contracted 17
         expanded 16 contracted 54
         expanded 39 contracted 17
         expanded 16 contracted 53
         expanded 39 contracted 17
         expanded 16 contracted 53
         expanded 39 contracted 17
         expanded 16 contracted 53
         expanded 39 contracted 17
         probability: 2.828, fragment count: 50779, expected: 40000
         
         expanded 0 contracted 4168
         expanded 7 contracted 505
         expanded 1 contracted 71
         expanded 7 contracted 1
         expanded 0 contracted 15
         expanded 7 contracted 0
         expanded 0 contracted 14
         expanded 7 contracted 0
         expanded 0 contracted 14
         expanded 7 contracted 0
         probability: 4.0, fragment count: 38109, expected: 40000
         expanded 713 contracted 0
         expanded 0 contracted 0
         expanded 0 contracted 14
         expanded 7 contracted 0
         expanded 0 contracted 14
         expanded 7 contracted 0
         expanded 0 contracted 14
         expanded 7 contracted 0
         expanded 0 contracted 14
         expanded 7 contracted 0
         expanded 0 contracted 14
         expanded 7 contracted 0
         probability: 2.828, fragment count: 50996, expected: 40000
         
         expanded 1 contracted 0
         expanded 8 contracted 0
         expanded 64 contracted 0
         expanded 512 contracted 0
         expanded 8 contracted 4088
         expanded 8 contracted 4112
         expanded 8 contracted 4136
         expanded 8 contracted 4152
         expanded 48 contracted 4152
         expanded 216 contracted 4152
         expanded 424 contracted 4152
         expanded 112 contracted 4168
         expanded 144 contracted 4168
         expanded 32 contracted 4168
         expanded 0 contracted 4168
         expanded 7 contracted 505
         expanded 1 contracted 71
         expanded 7 contracted 1
         probability: 4.0, fragment count: 38025, expected: 40000
         expanded 704 contracted 0
         expanded 16 contracted 14
         expanded 0 contracted 14
         expanded 7 contracted 0
         probability: 2.828, fragment count: 50912, expected: 40000
         
         expanded 1 contracted 0
         expanded 8 contracted 0
         expanded 64 contracted 0
         expanded 512 contracted 0
         expanded 8 contracted 4088
         expanded 8 contracted 4112
         expanded 8 contracted 4136
         expanded 8 contracted 4152
         expanded 48 contracted 4168
         expanded 216 contracted 4200
         expanded 424 contracted 4248
         expanded 112 contracted 4352
         expanded 144 contracted 4480
         expanded 32 contracted 4496
         expanded 0 contracted 4536
         expanded 143 contracted 564
         probability: 4.0, fragment count: 39376, expected: 40000
         expanded 758 contracted 504
         expanded 128 contracted 804
         expanded 0 contracted 876
         expanded 174 contracted 93
         probability: 2.828, fragment count: 48728, expected: 40000
         
         expanded 1 contracted 0
         expanded 8 contracted 0
         expanded 64 contracted 0
         exp
         anded 512 contracted 0
         expanded 8 contracted 4088
         expanded 8 contracted 4112
         expanded 8 contracted 4136
         expanded 8 contracted 4152
         expanded 48 contracted 4168
         expanded 216 contracted 4200
         expanded 424 contracted 4248
         expanded 112 contracted 4352
         expanded 144 contracted 4480
         expanded 32 contracted 4496
         expanded 0 contracted 4536
         expanded 143 contracted 564
         expanded 184 contracted 60
         expanded 185 contracted 2
         expanded 186 contracted 1
         expanded 186 contracted 1
         expanded 186 contracted 1
         expanded 186 contracted 0
         probability: 4.0, fragment count: 34973, expected: 40000
         expanded 757 contracted 0
         expanded 136 contracted 300
         expanded 0 contracted 372
         expanded 190 contracted 37
         expanded 225 contracted 5
         expanded 224 contracted 2
         expanded 224 contracted 1
         expanded 224 contracted 1
         expanded 224 contracted 0
         probability: 2.828, fragment count: 48112, expected: 40000
         */
        if thresholdMaskSum > 0 {
          expanded.append((UInt32(nodeID), thresholdMaskSum))
        } else if importanceMetric.sum() < threshold,
                  branchesMask == 0 {
          contracted.append(UInt32(nodeID))
        }
      }
      
      consoleMessage += "expanded \(expanded.count) contracted \(contracted.count)"
      consoleMessage += "\n"
      
      return (expanded, contracted)
    }
    
    func resizeOctreeNodes(
      expanded: [(UInt32, UInt8)],
      contracted: [UInt32]
    ) {
      let oldToNewMap = octree.resizeNodes(
        expanded: expanded, contracted: contracted)
      
      let count = octree.nodes.count
      let NAN = SIMD8<Float>(repeating: .nan)
      var newDensityIntegrals = Array(repeating: NAN, count: count)
      var newGradientIntegrals = Array(repeating: NAN, count: count)
      for oldNodeID in densityIntegrals.indices {
        let newNodeID = Int(oldToNewMap[oldNodeID])
        if newNodeID < UInt32.max {
          newDensityIntegrals[newNodeID] = densityIntegrals[oldNodeID]
          newGradientIntegrals[newNodeID] = gradientIntegrals[oldNodeID]
        }
      }
      densityIntegrals = newDensityIntegrals
      gradientIntegrals = newGradientIntegrals
    }
    
    // If we split at the maximum probability, the number of cells averages out
    // to roughly 4x the intended number. Therefore, we start at a looser
    // tolerance and repeat if the fragment count undershoots.
    let probabilityMultipliers: [Float] = [
      4, 2.828, 2, 1.414, 1
    ]
    for probabilityMultiplier in probabilityMultipliers {
      var iterationID = 0
      var enableContraction = false
      while true {
        let (expanded, contracted) = queryOctreeNodes(
          probabilityMultiplier: probabilityMultiplier)
        
        if expanded.count == 0 {
          if contracted.count == 0 {
            break
          } else {
            enableContraction = true
          }
        } else if contracted.count == 0 {
          if enableContraction {
            break
          }
        }
        
        if enableContraction {
          resizeOctreeNodes(expanded: [], contracted: contracted)
        } else {
          resizeOctreeNodes(expanded: expanded, contracted: [])
        }
        
        iterationID += 1
        if iterationID > 50 {
          var octreeFragmentCount = 0
          for node in octree.nodes {
            let leafMask = ~node.branchesMask
            octreeFragmentCount += leafMask.nonzeroBitCount
          }
          fatalError("""
            Wave function failed to converge after 50 iterations.
            Fragment Count: \(octreeFragmentCount)
            Probability Multiplier: \(probabilityMultiplier)
            """)
        }
      }
      
      var octreeFragmentCount = 0
      for node in octree.nodes {
        let leafMask = ~node.branchesMask
        octreeFragmentCount += leafMask.nonzeroBitCount
      }
      consoleMessage += "probability: \(probabilityMultiplier), fragment count: \(octreeFragmentCount), expected: \(fragmentCount)\n"
      
      if octreeFragmentCount >= fragmentCount {
        break
      } else if probabilityMultiplier == 1 {
        fatalError("Could not create octree with the specified fragment count.")
      }
    }
    print(consoleMessage)
    
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
  }
}
