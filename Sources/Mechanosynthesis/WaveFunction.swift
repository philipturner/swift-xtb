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
          // - Compute integrals with 8x higher accuracy.
          // - Report what fragment count is actually being used.
          // - Reduce the actual fragment count by 8x.
          // - Sample importance at 64x accuracy, fixing flerovium.
          // - Optimize the integral computation, ensuring it still
          //   produces the same results to within rounding error.
          
          // Lookup table for child nodes.
          var x = SIMD8<Float>(0, 1, 0, 1, 0, 1, 0, 1) * 0.5 - 0.25
          var y = SIMD8<Float>(0, 0, 1, 1, 0, 0, 1, 1) * 0.5 - 0.25
          var z = SIMD8<Float>(0, 0, 0, 0, 1, 1, 1, 1) * 0.5 - 0.25
          x = x * node.spacing + node.center.x
          y = y * node.spacing + node.center.y
          z = z * node.spacing + node.center.z
          
          for branchID in 0..<8 {
            let xyz = SIMD3(x[branchID], y[branchID], z[branchID])
            let Ψ = createCellValues(center: xyz, spacing: node.spacing / 2)
            let eighthNodeSpacing = node.spacing / 64
            let d3r = node.spacing * node.spacing * eighthNodeSpacing
            let _densityIntegral = Ψ * 1 * Ψ * d3r
            
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
            var _gradientIntegral = finalIntΔx + finalIntΔy + finalIntΔz
            _gradientIntegral *= eighthNodeSpacing
            
            densityIntegral[branchID] = _densityIntegral.sum()
            gradientIntegral[branchID] = _gradientIntegral.sum()
          }
          
          let ε: Float = .leastNormalMagnitude
          densityIntegral.replace(with: ε, where: densityIntegral .< ε)
          gradientIntegral.replace(with: ε, where: gradientIntegral .< ε)
          densityIntegrals[nodeID] = densityIntegral
          gradientIntegrals[nodeID] = gradientIntegral
        }
        
        if unsafeBitCast(node.branchesMask, to: UInt64.self) != .zero {
          var mask32 = SIMD8<UInt32>(
            truncatingIfNeeded: SIMD8(repeating: 1) &- node.branchesMask)
          mask32 &*= Float(1).bitPattern
          densityIntegral *= unsafeBitCast(mask32, to: SIMD8<Float>.self)
          gradientIntegral *= unsafeBitCast(mask32, to: SIMD8<Float>.self)
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
      expanded: [SIMD8<UInt8>]?,
      contracted: [Bool]?
    ) {
      let (globalDensity, globalSquareGradient) = evaluateGlobalIntegrals()
      let count = octree.nodes.count
      var expanded = [SIMD8<UInt8>](repeating: .zero, count: count)
      var contracted = [Bool](repeating: false, count: count)
      var expandedIsZero = true
      var contractedIsZero = true
      
      let threshold = probabilityMultiplier / Float(fragmentCount)
      let densityWeight = 0.5 / globalDensity
      let gradientWeight = 0.5 / globalSquareGradient
      
      for nodeID in octree.nodes.indices {
        let importanceMetric =
        densityIntegrals[nodeID] * densityWeight +
        gradientIntegrals[nodeID] * gradientWeight
        let branchesMask = octree.nodes[nodeID].branchesMask
        
        if isContracting {
          guard unsafeBitCast(branchesMask, to: UInt64.self) == .zero else {
            continue
          }
          
          // Soon, the threshold will be increased 8x, so this comparison will
          // represent the actual threshold.
          let importanceMax = importanceMetric.max()
          if importanceMax < threshold / 8 {
            contracted[nodeID] = true
            contractedIsZero = false
          }
        } else {
          var mask32 = SIMD8<UInt32>(repeating: .zero)
          mask32.replace(with: 1, where: importanceMetric .> threshold)
          var mask8 = SIMD8<UInt8>(truncatingIfNeeded: mask32)
          mask8 &= SIMD8(repeating: 1) &- branchesMask
          
          if unsafeBitCast(mask8, to: UInt64.self) != .zero {
            expanded[nodeID] = mask8
            expandedIsZero = false
          }
        }
      }
      
      return (
        expandedIsZero ? nil : expanded,
        contractedIsZero ? nil : contracted)
    }
    
    func resizeOctreeNodes(
      expanded: [SIMD8<UInt8>],
      contracted: [Bool]
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
          probabilityMultiplier: probabilityMultiplier,
          isContracting: enableContraction)
        
        if enableContraction {
          if contracted == nil {
            break
          } else {
            // Don't waste time if this round of contractions won't produce a
            // converged solution.
            var octreeFragmentCount = 0
            for node in octree.nodes {
              let mask64 = unsafeBitCast(node.branchesMask, to: UInt64.self)
              octreeFragmentCount += 8 - mask64.nonzeroBitCount
            }
            if octreeFragmentCount <= fragmentCount {
              break
            }
          }
        } else {
          if expanded == nil {
            enableContraction = true
          }
        }
        
        if enableContraction {
          if let contracted {
            let count = octree.nodes.count
            let expanded = [SIMD8<UInt8>](repeating: .zero, count: count)
            resizeOctreeNodes(expanded: expanded, contracted: contracted)
          }
        } else {
          let count = octree.nodes.count
          let contracted = [Bool](repeating: false, count: count)
          guard let expanded else {
            fatalError("This should never happen.")
          }
          resizeOctreeNodes(expanded: expanded, contracted: contracted)
        }
        
        iterationID += 1
        if iterationID > 50 {
          var octreeFragmentCount = 0
          for node in octree.nodes {
            let mask64 = unsafeBitCast(node.branchesMask, to: UInt64.self)
            octreeFragmentCount += 8 - mask64.nonzeroBitCount
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
        let mask64 = unsafeBitCast(node.branchesMask, to: UInt64.self)
        octreeFragmentCount += 8 - mask64.nonzeroBitCount
      }
      
      if octreeFragmentCount >= fragmentCount {
        break
      } else if probabilityMultiplier == 1 {
        fatalError("Could not create octree with the specified fragment count.")
      }
    }
  }
}
