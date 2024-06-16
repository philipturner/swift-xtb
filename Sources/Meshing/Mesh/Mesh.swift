//
//  Mesh.swift
//
//
//  Created by Philip Turner on 5/26/24.
//

public struct MeshDescriptor {
  // Required. The set of octrees to generate the topology from.
  public var octrees: [Octree] = []
  
  // Optional. The shift of each octree with respect to the origin.
  //
  // To initialize the ansatz easily and accurately, the octree should be
  // origin-centered and with minimal size exponent. This setting
  // allows different origin-centered meshes to be created for different
  // atoms, then fused after shifting to the actual position.
  //
  // Each shift must be a multiple of the coarse voxel spacing.
  public var positions: [SIMD3<Int32>]?
  
  // Required. The power-2 size of coarse voxels, in Bohr.
  public var sizeExponent: Int?
  
  public init() {
    
  }
}

public struct Mesh {
  // A uniform grid of coarse voxels.
  public var coarseVoxels: Grid<CoarseVoxel>
  
  // The separation between coarse voxels, in Bohr.
  public var spacing: Int
  
  public init(descriptor: MeshDescriptor) {
    // Check the correctness of the size exponent.
    guard let sizeExponent = descriptor.sizeExponent else {
      fatalError("Descriptor was incomplete.")
    }
    guard sizeExponent > 0 else {
      fatalError("Coarse voxel spacing must be at least 2 Bohr.")
    }
    
    // Set the coarse voxel spacing.
    spacing = 1 << sizeExponent
    Self.checkOctreeSizes(
      octrees: descriptor.octrees,
      spacing: spacing)
    
    // Detach the nodes from the octrees.
    let nodes = Self.detachOctreeNodes(
      octrees: descriptor.octrees,
      positions: descriptor.positions,
      spacing: spacing)
    
    // Create an empty grid with the smallest possible bounding box.
    let globalBoundingBox = Self.createGlobalBoundingBox(
      nodes: nodes,
      spacing: spacing)
    coarseVoxels = Self.createCoarseGrid(
      minimum: globalBoundingBox.minimum,
      maximum: globalBoundingBox.maximum)
    
    // Place the nodes into an array for each voxel.
    let coarseNodeArrays = Self.mapNodesToCoarseVoxels(
      nodes: nodes,
      spacing: spacing,
      coarseVoxels: coarseVoxels)
    
    // Next:
    // - Function that initializes all the bounding box accumulators that will
    //   be generated in parallel.
    //   - (h = 2, spacing = 4)
    //   - (h = 1, spacing = 2)
    //   - h = 1
    // - After reduction, it allocates memory for the respective arrays.
    //   - grid of SIMD8<Float>
    //   - grid of SIMD8<Float>
    //   - grid of FineVoxel
    // - Stores information about occupied cells in 'levels' by a mask of
    //   Float.NaN.
    for voxelID in coarseVoxels.cells.indices {
      // Create the coarse voxel.
      let nodes = coarseNodeArrays[voxelID]
      var coarseVoxel = Self.createCoarseVoxel(
        nodes: nodes,
        spacing: spacing)
      Self.fillCoarseLevels(
        nodes: nodes,
        coarseVoxel: &coarseVoxel)
      
      // Create the fine voxels within the coarse voxel.
      
      // Store in the grid of coarse voxels.
      coarseVoxels.cells[voxelID] = coarseVoxel
    }
  }
}
