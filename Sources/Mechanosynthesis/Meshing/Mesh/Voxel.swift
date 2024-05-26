//
//  Voxel.swift
//
//
//  Created by Philip Turner on 5/26/24.
//

// A section of real space larger than 1 cubic Bohr.
public struct CoarseVoxel {
  // A multigrid of coarse cells, arranged in Morton order in 2x2x2 chunks.
  //
  // The levels are ordered from coarsest to finest. Their offsets are relative
  // to the lower corner of this voxel.
  var coarseLevels: [Grid<SIMD8<Float>>]
}

// A section of real space equal to 1 cubic Bohr.
public struct FineVoxel {
  // A multigrid of fine cells, arranged in Morton order in 2x2x2 chunks.
  //
  // The levels are ordered from coarsest to finest. Their offsets are relative
  // to the lower corner of this voxel.
  var fineLevels: [Grid<SIMD8<Float>>]
}
