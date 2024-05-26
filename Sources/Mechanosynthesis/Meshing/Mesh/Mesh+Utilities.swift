//
//  Mesh+Sorting.swift
//
//
//  Created by Philip Turner on 5/26/24.
//

extension Mesh {
  // Creates an empty coarse grid from the descriptor.
  static func createCoarseGrid(
    descriptor: MeshDescriptor
  ) -> Grid<CoarseVoxel> {
    let spacing = 1 << descriptor.sizeExponent!
    
    // Find the bounding box of the coarse grid.
    var minimum: SIMD3<Int32> = .init(repeating: .max)
    var maximum: SIMD3<Int32> = .init(repeating: -.max)
    for octree in descriptor.octrees {
      // The highest octree level must span 2 * coarse voxel spacing,
      // so that either half of the octree is 1 * coarse voxel spacing.
      let expectedSpacing = Float(2 * spacing)
      guard octree.nodes[0].spacing >= expectedSpacing else {
        fatalError("Octree must be large enough to fit the coarse voxels.")
      }
      
      for node in octree.nodes {
        // There are some topologies where no nodes will meet this criterion.
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
    }
    guard all(minimum .< maximum) else {
      // This happens when there are no octrees.
      fatalError("Mesh bounds could not be established.")
    }
    
    // Create the coarse grid.
    return Self.createCoarseGrid(minimum: minimum, maximum: maximum)
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
  // Maps octree nodes to coarse voxels.
  //
  // Unsure how to proceed with fine voxels for now.
  func mapNodesToCoarseVoxels(octree: Octree) -> [SIMD2<UInt32>] {
    // Create an accumulator for each voxel, just for this octree.
    var voxelAccumulators = [UInt32](
      repeating: .zero, count: coarseVoxels.cells.count)
    
    // Iterate over the nodes.
    var nodesToCoarseVoxelsMap: [SIMD2<UInt32>] = []
    for node in octree.nodes {
      guard node.spacing <= Float(spacing) else {
        // Append something here, to preserve the index in the array of nodes.
        let mappedLocation = SIMD2<UInt32>(repeating: .max)
        nodesToCoarseVoxelsMap.append(mappedLocation)
        continue
      }
      
      // Locate this node within the grid of voxels.
      let coordinate = node.center / Float(spacing)
      let coordinateDelta = coordinate - SIMD3<Float>(coarseVoxels.offset)
      guard all(coordinateDelta .> 0),
            all(coordinateDelta .< SIMD3(coarseVoxels.dimensions)) else {
        fatalError("Voxel was out of bounds.")
      }
      var voxelID: UInt32 = .zero
      do {
        let deltaInt = SIMD3<UInt32>(coordinateDelta.rounded(.down))
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
    return nodesToCoarseVoxelsMap
  }
}
