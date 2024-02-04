//
//  MultiGrid.swift
//
//
//  Created by Philip Turner on 2/4/24.
//

import OpenCL

// Primitive data type used everywhere in the framework.
// - occupancy encoded into the 32 bits per cell
// - data sourced from other levels is marked with flag bits
// - 20 bits of mantissa (s1e8m19)
//   - 6 decimal places for wavefunction
//   - doesn't limit precision, because sums are in FP32, and every group of
//     ~32 elements is Kahan summed
//
// OpenCL representation:
//
// typedef struct {
//   uint address;
//   float spacing;
//   int3 origin;
//   uint3 size;
// } MultiGridLevel;
public struct MultiGridLevel {
  // RAM address (in 4-byte words)
  public var address: UInt32 = .zero
  
  // spacing (power-2 multiple of Bohr)
  public var spacing: Float = 1
  
  // origin (as even integer multiple of 'spacing' away from .zero)
  public var origin: SIMD3<Int32> = .zero
  
  // allocated size (as even integer)
  // - Not every voxel within this volume is iterated over. It is just an
  //   efficient way to manage the memory.
  public var size: SIMD3<UInt32> = .zero
  
  public init() {}
}
