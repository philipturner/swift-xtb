//
//  Mesh+Sorting.swift
//
//
//  Created by Philip Turner on 5/26/24.
//

extension Mesh {
  static func checkOctreeSizes(
    octrees: [Octree],
    spacing: Int
  ) {
    for octree in octrees {
      // The highest octree level must span 2 * coarse voxel spacing,
      // so that either half of the octree is 1 * coarse voxel spacing.
      let expectedSpacing = Float(2 * spacing)
      guard octree.nodes[0].spacing >= expectedSpacing else {
        fatalError("Octree must be large enough to fit the coarse voxels.")
      }
    }
  }
  
  // Transform the nodes into a common coordinate space, then detach them from
  // their parent octrees.
  static func detachOctreeNodes(
    octrees: [Octree],
    positions: [SIMD3<Int32>]?,
    spacing: Int
  ) -> [OctreeNode] {
    // Find the position of each octree in the global coordinate space.
    var octreePositions: [SIMD3<Int32>]
    if let positions {
      octreePositions = positions
    } else {
      octreePositions = Array(repeating: .zero, count: octrees.count)
    }
    guard octreePositions.count == octrees.count else {
      fatalError("Number of positions did not match number of octrees.")
    }
    for position in octreePositions {
      guard all(position % Int32(spacing) .== 0) else {
        fatalError("Octree position must be aligned to coarse voxel spacing.")
      }
    }
    
    // Collect up the nodes into a single array.
    var output: [OctreeNode] = []
    for octreeID in octrees.indices {
      let octree = octrees[octreeID]
      let octreePosition = octreePositions[octreeID]
      
      // Iterate over the nodes within this octree.
      for nodeID in octree.nodes.indices {
        // Shift the node's center.
        var node = octree.nodes[nodeID]
        var center = node.center
        center += SIMD3<Float>(octreePosition)
        node.centerAndSpacing = SIMD4(center, node.spacing)
        
        // Remove references to other nodes.
        node.parentIndex = .max
        node.branchesIndex = .max
        node.branchesMask.replace(
          with: SIMD8.zero, where: node.branchesMask .!= 255)
        output.append(node)
      }
    }
    return output
  }
  
  // Creates an empty coarse grid from the correctly shifted nodes.
  static func createGlobalBoundingBox(
    nodes: [OctreeNode],
    spacing: Int
  ) -> (
    minimum: SIMD3<Int32>,
    maximum: SIMD3<Int32>
  ) {
    // Find the bounding box of the coarse grid.
    var minimum: SIMD3<Int32> = .init(repeating: .max)
    var maximum: SIMD3<Int32> = .init(repeating: -.max)
    for node in nodes {
      // There are octree topologies where no nodes will meet this criterion.
      // For example, if the octree terminates with only a single node.
      //
      // Remember, the coarse voxel spacing is twice the resolution of the
      // actual data that is stored.
      guard node.spacing == Float(spacing) else {
        continue
      }
      
      // Find the bounding box of the node.
      let nodeMinimum = SIMD3<Int32>(node.center - node.spacing / 2)
      let nodeMaximum = SIMD3<Int32>(node.center + node.spacing / 2)
      guard all(nodeMinimum % Int32(spacing) .== 0),
            all(nodeMaximum % Int32(spacing) .== 0) else {
        fatalError("Node bounds were not aligned to spacing.")
      }
      
      // Find the nearest integer multiple of the voxel spacing.
      let alignedMinimum = nodeMinimum / Int32(spacing)
      let alignedMaximum = nodeMaximum / Int32(spacing)
      minimum.replace(with: alignedMinimum, where: alignedMinimum .< minimum)
      maximum.replace(with: alignedMaximum, where: alignedMaximum .> maximum)
    }
    guard all(minimum .< maximum) else {
      // This happens when there are no octrees.
      fatalError("Mesh bounds could not be established.")
    }
    
    // Return the bounding box.
    return (minimum, maximum)
  }
  
  // Creates an empty coarse grid from the bounding box.
  static func createCoarseGrid(
    minimum: SIMD3<Int32>, maximum: SIMD3<Int32>
  ) -> Grid<CoarseVoxel> {
    // Create an empty grid of fine voxels.
    var fineGridDesc = GridDescriptor<FineVoxel>()
    fineGridDesc.dimensions = .zero
    fineGridDesc.emptyElement = FineVoxel()
    let emptyFineGrid = Grid(descriptor: fineGridDesc)
    
    // Create an empty coarse voxel.
    var coarseVoxelDesc = CoarseVoxelDescriptor()
    coarseVoxelDesc.fineVoxels = emptyFineGrid
    let emptyCoarseVoxel = CoarseVoxel(descriptor: coarseVoxelDesc)
    
    // Create a grid of coarse voxels.
    var coarseGridDesc = GridDescriptor<CoarseVoxel>()
    coarseGridDesc.offset = minimum
    coarseGridDesc.dimensions = SIMD3(truncatingIfNeeded: maximum &- minimum)
    coarseGridDesc.emptyElement = emptyCoarseVoxel
    return Grid(descriptor: coarseGridDesc)
  }
}

extension Mesh {
  // Map octree nodes to coarse voxels.
  static func mapNodesToCoarseVoxels(
    nodes: [OctreeNode],
    spacing: Int,
    coarseVoxels: Grid<CoarseVoxel>
  ) -> [[OctreeNode]] {
    // Create an accumulator for each voxel.
    var voxelAccumulators = [UInt32](
      repeating: .zero, count: coarseVoxels.cells.count)
    
    // Find the indices for each node.
    var nodesToCoarseVoxelsMap: [SIMD2<UInt32>] = []
    for node in nodes {
      guard node.spacing <= Float(spacing) else {
        // Append something here, to preserve the index in the array of nodes.
        let mappedLocation = SIMD2<UInt32>(repeating: .max)
        nodesToCoarseVoxelsMap.append(mappedLocation)
        continue
      }
      
      // Locate this node within the grid of voxels.
      var voxelID: UInt32 = .zero
      do {
        var voxelCoordinate = node.center / Float(spacing)
        voxelCoordinate -= SIMD3<Float>(coarseVoxels.offset)
        guard all(voxelCoordinate .> 0),
              all(voxelCoordinate .< SIMD3(coarseVoxels.dimensions)) else {
          fatalError("Voxel was out of bounds.")
        }
        voxelCoordinate.round(.down)
        
        let deltaInt = SIMD3<UInt32>(voxelCoordinate)
        let dimensions = coarseVoxels.dimensions
        voxelID += deltaInt[0]
        voxelID += deltaInt[1] * dimensions[0]
        voxelID += deltaInt[2] * dimensions[0] * dimensions[1]
      }
      
      // Accumulate the number of nodes mapped to this voxel.
      let slotID = voxelAccumulators[Int(voxelID)]
      voxelAccumulators[Int(voxelID)] = slotID + 1
      let mappedLocation = SIMD2(voxelID, slotID)
      nodesToCoarseVoxelsMap.append(mappedLocation)
    }
    
    // Create an array for each voxel.
    var voxelArrays: [[OctreeNode]] = []
    for voxelID in coarseVoxels.cells.indices {
      let nodeCount = voxelAccumulators[voxelID]
      let emptyNode = OctreeNode(
        centerAndSpacing: SIMD4(repeating: .nan),
        parentIndex: UInt32.max,
        branchesIndex: UInt32.max,
        branchesMask: SIMD8(repeating: UInt8.max))
      let emptyArray = Array(repeating: emptyNode, count: Int(nodeCount))
      voxelArrays.append(emptyArray)
    }
    
    // Move each node to its corresponding voxel.
    for nodeID in nodes.indices {
      let map = nodesToCoarseVoxelsMap[nodeID]
      let voxelID = map[0]
      
      guard voxelID < .max else {
        // This node is too large to fit into the mesh.
        continue
      }
      
      // Fetch the node for modification and translocation.
      var node = nodes[nodeID]
      do {
        // Translate the node by the grid offset (in atomic units).
        var center = node.center
        center -= SIMD3<Float>(coarseVoxels.offset) * Float(spacing)
        
        // Find the cell address (in multiples of coarse voxel spacing).
        var voxelCoordinate = center / Float(spacing)
        voxelCoordinate.round(.down)
        
        // Scale the cell address back into atomic units, then subtract from
        // the center.
        center -= voxelCoordinate * Float(spacing)
        node.centerAndSpacing = SIMD4(center, node.spacing)
        guard all(node.center .> 0),
              all(node.center .< Float(spacing)) else {
          fatalError("Node was not fully contained inside coarse voxel.")
        }
      }
      
      // Write the node to the appropriate array slot.
      let slotID = map[1]
      voxelArrays[Int(voxelID)][Int(slotID)] = node
    }
    
    // Check that all slots have been filled.
    // - This will likely be a temporary check, removed after comprehensive
    //   testing.
    // - For now, it makes it easier to be sure that the meshing program is
    //   running correctly.
    for i in voxelArrays.indices {
      for j in voxelArrays[i].indices {
        let node = voxelArrays[i][j]
        if node.center.x.isNaN {
          fatalError("The nodes were not written correctly.")
        }
      }
    }
    
    return voxelArrays
  }
  
  // Allocate the coarse levels (1 ≤ h < coarse voxel size) for the
  // nodes. Use these allocations to create a new coarse voxel with the correct
  // grid sizes.
  static func createCoarseVoxel(
    nodes: [OctreeNode],
    spacing: Int
  ) -> CoarseVoxel {
    // Decode the size exponent from the mesh spacing.
    //
    // sizeExponent = 0
    // error
    //
    // sizeExponent = 1
    // array of FineVoxel / spacing = 1, count ≤ 8
    // array of h = 1     / spacing = 2, count ≤ 1
    //
    // sizeExponent = 2
    // array of FineVoxel / spacing = 1, count ≤ 64
    // array of h = 1     / spacing = 2, count ≤ 8
    // array of h = 2     / spacing = 4, count ≤ 1
    //
    // sizeExponent = 3
    // array of FineVoxel / spacing = 1, count ≤ 512
    // array of h = 1     / spacing = 2, count ≤ 64
    // array of h = 2     / spacing = 4, count ≤ 8
    // array of h = 4     / spacing = 8, count ≤ 1
    let sizeExponent = spacing.trailingZeroBitCount
    
    // Allocate bounding box accumulators for each level tracked.
    var levelMinima = [SIMD3<Int32>](
      repeating: SIMD3(repeating: .max), count: 1 + sizeExponent)
    var levelMaxima = [SIMD3<Int32>](
      repeating: SIMD3(repeating: -.max), count: 1 + sizeExponent)
    
    // Iterate over all the voxels with an integer spacing.
    for node in nodes {
      guard node.spacing >= 1 else {
        continue
      }
      
      // Find the bounding box of the node.
      let nodeMinimum = SIMD3<Int32>(node.center - node.spacing / 2)
      let nodeMaximum = SIMD3<Int32>(node.center + node.spacing / 2)
      let levelID = Int32(node.spacing).trailingZeroBitCount
      
      // Update the level's bounding box.
      var levelMinimum = levelMinima[levelID]
      var levelMaximum = levelMaxima[levelID]
      levelMinimum.replace(
        with: nodeMinimum, where: nodeMinimum .< levelMinimum)
      levelMaximum.replace(
        with: nodeMaximum, where: nodeMaximum .> levelMaximum)
      levelMinima[levelID] = levelMinimum
      levelMaxima[levelID] = levelMaximum
    }
    
    // Iterate over the levels.
    for levelID in 0..<(1 + sizeExponent) {
      let minimum = levelMinima[levelID]
      let maximum = levelMaxima[levelID]
      
      // If the level is empty , set the bounding box to zero.
      if any(minimum .> maximum) {
        levelMinima[levelID] = .zero
        levelMaxima[levelID] = .zero
      }
    }
    
    // Create a grid of fine voxels.
    var fineGridDesc = GridDescriptor<FineVoxel>()
    fineGridDesc.offset = levelMinima[0]
    fineGridDesc.dimensions = SIMD3(
      truncatingIfNeeded: levelMaxima[0] &- levelMinima[0])
    fineGridDesc.emptyElement = FineVoxel()
    
    // Create a coarse voxel.
    var coarseVoxelDesc = CoarseVoxelDescriptor()
    coarseVoxelDesc.fineVoxels = Grid(descriptor: fineGridDesc)
    var coarseVoxel = CoarseVoxel(descriptor: coarseVoxelDesc)
    
    // Iterate over the levels stored as floating-point numbers.
    for levelID in 1..<(1 + sizeExponent) {
      let minimum = levelMinima[levelID]
      let maximum = levelMaxima[levelID]
      
      // Create a grid of cells.
      var levelDesc = GridDescriptor<SIMD8<Float>>()
      levelDesc.offset = minimum
      levelDesc.dimensions = SIMD3(truncatingIfNeeded: maximum &- minimum)
      levelDesc.emptyElement = SIMD8(repeating: .nan)
      let level = Grid(descriptor: levelDesc)
      
      // Append to the levels in ascending order of voxel size.
      coarseVoxel.coarseLevels.append(level)
    }
    
    return coarseVoxel
  }
  
  // Fill the data for the coarse levels within a voxel.
  static func fillCoarseLevels(
    nodes: [OctreeNode],
    coarseVoxel: inout CoarseVoxel
  ) {
    
  }
}
