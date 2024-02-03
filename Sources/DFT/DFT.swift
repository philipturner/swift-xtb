
// Where to start?
// - Pick anything that's a reasonable name for a wavefunction data structure.
// - Lay out the sparse internal representation that will eventually be ported
//   to GPU.
// - Add some unit tests that query the expectation values of the different
//   hydrogen orbitals.

import OpenCL

struct MultiGridLevelDescriptor {
  var origin: SIMD3<Int> = .zero
  
  var spacing: Float = 1
  
  var size: SIMD3<Int> = .zero
  
  // perhaps offset within the larger OpenCL buffer allocation, which removes
  // the need to create a separate CLBuffer for every allocation (cheaper)
  
  init() {
    
  }
}

// Primitive data type used everywhere in the framework. Reallocated quite
// frequently.
// - occupancy encoded into the 32 bits per cell
// - data sourced from other levels is marked with flag bits
// - 20 bits of mantissa (6 decimal places)
struct MultiGridLevel {
  // spacing (power of 2 Bohr)
  // origin (as even integer multiple of 'spacing' away from .zero)
  // size (as even integer)
  // OpenCL buffer of 'float'
  //
  // OpenCL representation:
  //
  // typedef struct {
  //   int3 origin;
  //   float spacing;
  //   uint3 size;
  //   global uint* contents; // may required a preceding GPU kernel to force
  //                             the GPU address into the structure?
  //                          // perhaps create a massive 4 GB buffer and
  //                             sub-allocate from that on the CPU, making
  //                             pointers 32-bit
  //                          // TODO: can we fetch the GPU address via a GPU
  //                             kernel, read on the CPU, and insert into data
  //                             structures for a subsequent GPU kernel?
  // } MultiGridLevel;
  init(descriptor: MultiGridLevelDescriptor) {
    
  }
}
