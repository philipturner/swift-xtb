//
//  HydrogenicOrbital.swift
//
//
//  Created by Philip Turner on 2/13/24.
//

// A configuration for a hydrogenic orbital.
struct HydrogenicOrbitalDescriptor {
  // The functional form of the initial guess.
  var basisFunction: BasisFunction?
  
  // The minimum number of fragments to split each electron into.
  var fragmentCount: Int?
  
  // The nuclear position for the initial guess.
  var nucleusPosition: SIMD3<Float>?
  
  // The power-2 size of the coarsest octree level.
  var sizeExponent: Int?
}

/// A hydrogenic orbital, with a topology that divides real space into regions
/// of high charge density.
public struct HydrogenicOrbital {
  /// The octree that stores the structure of the wavefunction.
  public var octree: Octree
  
  /// The basis function for the orbital.
  public var basisFunction: BasisFunction
  
  init(descriptor: HydrogenicOrbitalDescriptor) {
    guard let basisFunction = descriptor.basisFunction,
          let fragmentCount = descriptor.fragmentCount,
          let nucleusPosition = descriptor.nucleusPosition,
          let sizeExponent = descriptor.sizeExponent else {
      fatalError("Descriptor was incomplete.")
    }
    self.basisFunction = basisFunction
    
    var octreeDesc = OctreeDescriptor()
    octreeDesc.sizeExponent = sizeExponent
    self.octree = Octree(descriptor: octreeDesc)
    
    // Cache the integrals so they don't need to be recomputed.
    var densityIntegrals: [SIMD8<Float>] = []
    var gradientIntegrals: [SIMD8<Float>] = []
    densityIntegrals.append(SIMD8(repeating: .nan))
    gradientIntegrals.append(SIMD8(repeating: .nan))
    
    @_transparent
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
      
      // Omit the angular part to ensure the mesh is rotationally invariant.
      //
      // In addition, omit the scaling by 1 / (4 * Float.pi).squareRoot().
      // This optimization slightly reduces the compute cost, without
      // changing the relative weighting of the cells.
      
      // NOTE: There is an opportunity to reduce the compute cost of ansatz
      // generation. Orbitals within the same sub-shell have the same mesh.
      let r = (x * x + y * y + z * z).squareRoot()
      return basisFunction.radialPart(r: r)
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
          // Lookup table for child nodes.
          var lx = SIMD8<Float>(0, 1, 0, 1, 0, 1, 0, 1) * 0.5 - 0.25
          var ly = SIMD8<Float>(0, 0, 1, 1, 0, 0, 1, 1) * 0.5 - 0.25
          var lz = SIMD8<Float>(0, 0, 0, 0, 1, 1, 1, 1) * 0.5 - 0.25
          lx = lx * node.spacing + node.center.x
          ly = ly * node.spacing + node.center.y
          lz = lz * node.spacing + node.center.z
          
          for branchID in 0..<8 {
            var x = SIMD8<Float>(0, 1, 0, 1, 0, 1, 0, 1) * 0.5 - 0.25
            var y = SIMD8<Float>(0, 0, 1, 1, 0, 0, 1, 1) * 0.5 - 0.25
            var z = SIMD8<Float>(0, 0, 0, 0, 1, 1, 1, 1) * 0.5 - 0.25
            x = x * node.spacing / 2 + lx[branchID]
            y = y * node.spacing / 2 + ly[branchID]
            z = z * node.spacing / 2 + lz[branchID]
            
            // Sample the wavefunction at (8x)->(64x) resolution to produce
            // higher-quality integrals. This reduces the number of fragments
            // needed for core electrons of high-Z elements.
            var _densityIntegral: SIMD8<Float> = .zero
            var _gradientIntegral: SIMD8<Float> = .zero
            for subBranchID in 0..<8 {
              let xyz = SIMD3(x[subBranchID], y[subBranchID], z[subBranchID])
              let Ψ = createCellValues(center: xyz, spacing: node.spacing / 4)
              let d3r = node.spacing * node.spacing * node.spacing / 512
              let __densityIntegral = Ψ * 1 * Ψ * d3r
              _densityIntegral[subBranchID] = __densityIntegral.sum()
              
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
              var __gradientIntegral: SIMD4<Float> = .zero
              __gradientIntegral = Δx * Δx + Δy * Δy + Δz * Δz
              _gradientIntegral[subBranchID] = __gradientIntegral.sum()
            }
            densityIntegral[branchID] = _densityIntegral.sum()
            gradientIntegral[branchID] = _gradientIntegral.sum()
          }
          gradientIntegral *= 2 * node.spacing / 512
          
          let ε: Float = .leastNormalMagnitude
          densityIntegral.replace(with: ε, where: densityIntegral .< ε)
          gradientIntegral.replace(with: ε, where: gradientIntegral .< ε)
          densityIntegrals[nodeID] = densityIntegral
          gradientIntegrals[nodeID] = gradientIntegral
        }
        
        if unsafeBitCast(node.branchesMask, to: UInt64.self) != .max {
          var mask32 = SIMD8<UInt32>(
            truncatingIfNeeded: node.branchesMask &>> 7)
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
          guard unsafeBitCast(branchesMask, to: UInt64.self) == .max else {
            continue
          }
          
          // 1x1x1 cell coordinating compute work
          // 2x2x2 cell < threshold / 8
          // 4x4x4 cell < threshold / 64
          let importanceMax = importanceMetric.max()
          if importanceMax < threshold / 8 {
            contracted[nodeID] = true
            contractedIsZero = false
          }
        } else {
          // 1x1x1 cell coordinating compute work
          // 2x2x2 cell > threshold
          // 4x4x4 cell > threshold / 8
          var mask32 = SIMD8<UInt32>(repeating: .zero)
          mask32.replace(with: 1, where: importanceMetric .> threshold)
          var mask8 = SIMD8<UInt8>(truncatingIfNeeded: mask32)
          mask8 &= branchesMask &>> 7
          
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
              let mask8 = node.branchesMask & SIMD8(repeating: 128)
              let mask64 = unsafeBitCast(mask8, to: UInt64.self)
              octreeFragmentCount += mask64.nonzeroBitCount
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
            let mask8 = node.branchesMask & SIMD8(repeating: 128)
            let mask64 = unsafeBitCast(mask8, to: UInt64.self)
            octreeFragmentCount += mask64.nonzeroBitCount
          }
          fatalError("""
            Meshing failed to converge after 50 iterations.
            Fragment Count: \(octreeFragmentCount)
            Probability Multiplier: \(probabilityMultiplier)
            """)
        }
      }
      
      var octreeFragmentCount = 0
      for node in octree.nodes {
        let mask8 = node.branchesMask & SIMD8(repeating: 128)
        let mask64 = unsafeBitCast(mask8, to: UInt64.self)
        octreeFragmentCount += mask64.nonzeroBitCount
      }
      
      if octreeFragmentCount >= fragmentCount {
        break
      } else if probabilityMultiplier == 1 {
        fatalError("Could not create octree with the specified fragment count.")
      }
    }
  }
}
