//
//  Voxel.swift
//
//
//  Created by Philip Turner on 5/26/24.
//

public struct CoarseVoxelDescriptor {
  // Required. A uniform grid of fine voxels.
  //
  // The grid's offset is relative to the lower corner of this voxel.
  public var fineVoxels: Grid<FineVoxel>?
}

// A large region of real space.
public struct CoarseVoxel {
  // A multigrid of 2x2x2 chunks of cells.
  //
  // The levels are ordered from finest to coarsest. Their offsets are relative
  // to the lower corner of this voxel.
  public var coarseLevels: [Grid<SIMD8<Float>>] = []
  
  // A uniform grid of fine voxels.
  public var fineVoxels: Grid<FineVoxel>
  
  public init(descriptor: CoarseVoxelDescriptor) {
    guard let fineVoxels = descriptor.fineVoxels else {
      fatalError("Descriptor was incomplete.")
    }
    self.fineVoxels = fineVoxels
  }
}

// A region of real space spanning 1 cubic Bohr.
public struct FineVoxel {
  // A multigrid of 2x2x2 chunks of cells.
  //
  // The levels are ordered from coarsest to finest. Their offsets are relative
  // to the lower corner of this voxel.
  public var fineLevels: [Grid<SIMD8<Float>>] = []
}
