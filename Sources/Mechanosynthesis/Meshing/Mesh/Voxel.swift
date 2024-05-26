//
//  Voxel.swift
//
//
//  Created by Philip Turner on 5/26/24.
//

public struct CoarseVoxelDescriptor {
  // Required. A uniform grid of fine voxels.
  //
  // The grid's offset is relative to the lower corner of this voxel. The
  // bounding box must fall within [0, 2^sizeExponent).
  public var fineVoxels: Grid<FineVoxel>?
  
  // Required. The power-2 size of the voxel, in Bohr.
  public var sizeExponent: Int?
}

// A large region of real space.
public struct CoarseVoxel {
  // A multigrid of 2x2x2 chunks of cells.
  //
  // The levels are ordered from coarsest to finest. Their offsets are relative
  // to the lower corner of this voxel.
  public var coarseLevels: [Grid<SIMD8<Float>>] = []
  
  // The size of the spanned region of real space.
  public var spacing: Int
}

// A region of real space spanning 1 cubic Bohr.
public struct FineVoxel {
  // A multigrid of 2x2x2 chunks of cells.
  //
  // The levels are ordered from coarsest to finest. Their offsets are relative
  // to the lower corner of this voxel.
  public var fineLevels: [Grid<SIMD8<Float>>] = []
}
