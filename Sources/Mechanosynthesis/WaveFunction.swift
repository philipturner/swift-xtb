//
//  WaveFunction.swift
//
//
//  Created by Philip Turner on 2/13/24.
//

// TODO: Remove this import!
import Foundation

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
          let eighthNodeSpacing = node.spacing / 8
          let d3r = node.spacing * node.spacing * eighthNodeSpacing
          densityIntegral = Ψ * 1 * Ψ * d3r
          
          let casted = unsafeBitCast(Ψ, to: SIMD4<UInt64>.self)
          let lowX = Ψ.evenHalf
          let highX = Ψ.oddHalf
          let lowY = unsafeBitCast(casted.evenHalf, to: SIMD4<Float>.self)
          let highY = unsafeBitCast(casted.oddHalf, to: SIMD4<Float>.self)
          let lowZ = Ψ.lowHalf
          let highZ = Ψ.highHalf
          
          let Δx = highX - lowX
          let Δy = highY - lowY
          let Δz = highZ - lowZ
          let intΔx = Δx * Δx
          let intΔy = Δy * Δy
          let intΔz = Δz * Δz
          
          let finalIntΔx = SIMD8<Float>(
            intΔx[0], intΔx[0], intΔx[1], intΔx[1],
            intΔx[2], intΔx[2], intΔx[3], intΔx[3])
          let castedIntΔy = unsafeBitCast(intΔy, to: SIMD2<UInt64>.self)
          let castedIntΔy2 = SIMD4(castedIntΔy[0], castedIntΔy[0],
                                   castedIntΔy[1], castedIntΔy[1])
          let finalIntΔy = unsafeBitCast(castedIntΔy2, to: SIMD8<Float>.self)
          let finalIntΔz = SIMD8<Float>(lowHalf: intΔz, highHalf: intΔz)
          gradientIntegral = finalIntΔx + finalIntΔy + finalIntΔz
          gradientIntegral *= eighthNodeSpacing
          
          let ε: Float = .leastNormalMagnitude
          densityIntegral.replace(with: ε, where: densityIntegral .< ε)
          gradientIntegral.replace(with: ε, where: gradientIntegral .< ε)
          densityIntegrals[nodeID] = densityIntegral
          gradientIntegrals[nodeID] = gradientIntegral
        }
        
        if node.branchesMask != 0 {
          let masterMask = SIMD8<UInt8>(
            1 << 0, 1 << 1, 1 << 2, 1 << 3,
            1 << 4, 1 << 5, 1 << 6, 1 << 7)
          let mask8 = masterMask & node.branchesMask
          let mask32 = SIMD8<UInt32>(truncatingIfNeeded: mask8) .!= 0
          densityIntegral.replace(with: 0, where: mask32)
          gradientIntegral.replace(with: 0, where: mask32)
        }
        density += Double(densityIntegral.sum())
        squareGradient += Double(gradientIntegral.sum())
      }
      return (Float(density), Float(squareGradient))
    }
    
    // Perform one iteration of octree resizing.
    func queryOctreeNodes(
      probabilityMultiplier: Float,
      isContracting: Bool
    ) -> (
      expanded: [(UInt32, UInt8)],
      contracted: [UInt32]
    ) {
      let (globalDensity, globalSquareGradient) = evaluateGlobalIntegrals()
      var expanded: [(UInt32, UInt8)] = []
      var contracted: [UInt32] = []
      
      let threshold = probabilityMultiplier / Float(fragmentCount)
      let densityWeight = 0.5 / globalDensity
      let gradientWeight = 0.5 / globalSquareGradient
      
      for nodeID in octree.nodes.indices {
        let importanceMetric =
        densityIntegrals[nodeID] * densityWeight +
        gradientIntegrals[nodeID] * gradientWeight
        let branchesMask = octree.nodes[nodeID].branchesMask
        
        if isContracting {
          guard branchesMask == 0 else {
            continue
          }
          
          let importanceSum = importanceMetric.sum()
          if importanceSum < threshold {
            contracted.append(UInt32(nodeID))
          }
        } else {
          let masterMask = SIMD8<UInt8>(
            1 << 0, 1 << 1, 1 << 2, 1 << 3,
            1 << 4, 1 << 5, 1 << 6, 1 << 7)
          let masterMask32 = SIMD8<UInt32>(truncatingIfNeeded: masterMask)
          
          var thresholdMask32 = SIMD8<UInt32>(repeating: .zero)
          thresholdMask32.replace(
            with: masterMask32, where: importanceMetric .> threshold)
          var thresholdMask8 = SIMD8<UInt8>(truncatingIfNeeded: thresholdMask32)
          if branchesMask != 0 {
            thresholdMask8.replace(
              with: 0, where: (masterMask & branchesMask) .!= 0)
          }
          
          let thresholdMaskSum = thresholdMask8.wrappedSum()
          if thresholdMaskSum > 0 {
            expanded.append((UInt32(nodeID), thresholdMaskSum))
          }
        }
      }
      
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
    var timeElapsed0: Double = 0
    var timeElapsed1: Double = 0
    for probabilityMultiplier in probabilityMultipliers {
      var iterationID = 0
      var enableContraction = false
      while true {
        let checkpoint0 = Date()
        let (expanded, contracted) = queryOctreeNodes(
          probabilityMultiplier: probabilityMultiplier,
          isContracting: enableContraction)
        let checkpoint1 = Date()
        
        if enableContraction {
          if contracted.count == 0 {
            break
          }
        } else {
          if expanded.count == 0 {
            enableContraction = true
          }
        }
        
        if enableContraction {
          if contracted.count > 0 {
            resizeOctreeNodes(expanded: [], contracted: contracted)
          }
        } else {
          resizeOctreeNodes(expanded: expanded, contracted: [])
        }
        let checkpoint2 = Date()
        if octree.nodes.count > 10000 {
          timeElapsed0 += checkpoint1.timeIntervalSince(checkpoint0)
          timeElapsed1 += checkpoint2.timeIntervalSince(checkpoint1)
          let proportion0 = timeElapsed0 / (timeElapsed0 + timeElapsed1)
          let proportion1 = timeElapsed1 / (timeElapsed0 + timeElapsed1)
          print("iteration \(iterationID): \(proportion0), \(proportion1)")
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
      
      if octreeFragmentCount >= fragmentCount {
        break
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
  }
}
