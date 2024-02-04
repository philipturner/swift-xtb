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

// First goal: initialize wavefunctions in variable-resolution orbitals.
// - Method to initialize atomic orbitals on the GPU, as an ansatz to the
//   ground-state wavefunction.
// - Can be debugged before the multigrid data structure / level transfer
//   protocols are finalized.
// - Establish a procedure where the GPU recursively dispatches more threads in
//   subsequent kernels, entirely autonomous from the CPU. Use a fixed amount of
//   GPU threads (OpenMM heuristics), but atomics at the threadgroup scope
//   to achieve perfect saturation.
//   - user specifies the number of GPU threads upon SCF initialization
//     - better control for running multiple simulations simultaneously
//     - GPU kernel has no knowledge of the thread count; no need to encode it
//       into the kernel or embed it into the RAM
//   - 128-threadgroup granularity with 8-thread local memory seems reasonable
//   - the local memory protocol can easily be upgraded to subgroup shuffles
//   - single global atomic counter that's accessed serially
//   - exception: block summation happens at SIMD execution width granularity
// - Proactively encode GPU commands for a fixed number of multigrid levels,
// - Have the GPU automatically generate a multigrid for an atomic orbital, or
//   bail out if it runs out of memory.

struct MultiGrid {
  var levels: UnsafeMutablePointer<MultiGridLevel> // pointer to GPU RAM
}
