
// Where to start?
// - Pick anything that's a reasonable name for a wavefunction data structure.
// - Lay out the sparse internal representation that will eventually be ported
//   to GPU.
// - Add some unit tests that query the expectation values of the different
//   hydrogen orbitals.

import OpenCL

// TODO: Rewrite the multigrid API. The descriptor now stores the same data
// as the actual object. Try to avoid the descriptor paradigm in this library.
// You can initialize the 'SelfConsistentField' in a declarative manner, only
// requiring an OpenCL context in the initializer. The RAM can be initialized
// to 8 bytes. This API design is more expressive.

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
  //                             - expand to 16 GB and user-specified allocation
  //                               size by increasing computer word size from
  //                               8 bits to 32 bits
  //                             - the first argument of every GPU kernel is
  //                               global uint* RAM
  //                          // TODO: can we fetch the GPU address via a GPU
  //                             kernel, read on the CPU, and insert into data
  //                             structures for a subsequent GPU kernel?
  // } MultiGridLevel;
  init(descriptor: MultiGridLevelDescriptor) {
    
  }
}

public struct RAM {
  public var buffer: CLBuffer
  public var cpuAddress: UnsafeMutableRawPointer
  public var gpuAddress: UInt64
  
  public init(context: CLContext, size: Int) {
    let buffer = CLBuffer(context: context, flags: .readWrite, size: size)
    guard let buffer else {
      fatalError("Could not initialize OpenCL buffer.")
    }
    self.buffer = buffer
    
    // Create a temporary command queue.
    guard let device = context.devices?.first,
          let queue = CLCommandQueue(context: context, device: device) else {
      fatalError("Could not create queue.")
    }
    self.cpuAddress = try! queue.enqueueMap(
      buffer, flags: [.read, .write], offset: 0, size: size)
    
    // Compile the kernel for writing the GPU address.
    let source = """
    typedef struct {
      global uint* RAM;
    } PointerCapsule;
    
    // Store the GPU address in the first few bytes of RAM.
    kernel void copyAddress(global uint* RAM) {
      PointerCapsule capsule;
      capsule.RAM = RAM;
      ((global PointerCapsule*)RAM)[0] = capsule;
    }
    """
    guard let program = CLProgram(context: context, source: source) else {
      fatalError("Could not create program.")
    }
    do {
      try program.build()
    } catch {
      let log = program.buildLog(device: device)
      fatalError("Build error: \(log ?? "n/a")")
    }
    guard let kernels = program.createKernels(), kernels.count == 1 else {
      fatalError("Could not create kernel.")
    }
    
    // Encode the kernel.
    let kernel = kernels[0]
    try! kernel.setArgument(buffer, index: 0)
    try! queue.enqueueKernel(kernel, globalSize: CLNDRange(width: 1))
    try! queue.finish()
    
    let gpuAddressPointer = cpuAddress.assumingMemoryBound(to: UInt64.self)
    gpuAddress = gpuAddressPointer.pointee
  }
}

// struct SelfConsistentField - like a context, coordinates everything
