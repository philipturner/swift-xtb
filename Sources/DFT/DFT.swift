
// Where to start?
// - Pick anything that's a reasonable name for a wavefunction data structure.
// - Lay out the sparse internal representation that will eventually be ported
//   to GPU.
// - Add some unit tests that query the expectation values of the different
//   hydrogen orbitals.

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

public struct RAM {
  public var buffer: CLBuffer
  public var cpuAddress: UnsafeMutableRawPointer
  public var gpuAddress: UInt64
  
  public init(queue: CLCommandQueue, size: Int) {
    // 16 GB is the maximum around of 4-byte words addressable with 32 bits.
    guard size >= 0, size <= 16 * 1024 * 1024 * 1024 else {
      fatalError("Size was invalid.")
    }
    let buffer = CLBuffer(
      context: queue.context!, flags: .readWrite, size: size)
    guard let buffer else {
      fatalError("Could not initialize OpenCL buffer.")
    }
    self.buffer = buffer
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
    guard let program = CLProgram(
      context: queue.context!, source: source) else {
      fatalError("Could not create program.")
    }
    do {
      try program.build()
    } catch {
      let log = program.buildLog(device: queue.device!)
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

public class SelfConsistentField {
  @usableFromInline
  var _queue: CLCommandQueue?
  
  @usableFromInline
  var _ram: RAM?
  
  public init() {
    
  }
}

extension SelfConsistentField {
  public var queue: CLCommandQueue {
    @_transparent
    get {
      return _queue!
    }
    @_transparent
    set {
      _queue = newValue
    }
  }
  
  public var ram: RAM {
    @_transparent
    get {
      return _ram!
    }
    @_transparent
    set {
      _ram = newValue
    }
  }
}
