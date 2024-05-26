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
    guard let sizeExponent = descriptor.sizeExponent else {
      fatalError("Descriptor was incomplete.")
    }
    guard sizeExponent > 0 else {
      fatalError("Coarse voxel spacing must be at least 2 Bohr.")
    }
    spacing = 1 << sizeExponent
    
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
    coarseVoxels = Grid(descriptor: coarseGridDesc)
  }
}
