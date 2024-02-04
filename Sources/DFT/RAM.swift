//
//  RAM.swift
//
//
//  Created by Philip Turner on 2/4/24.
//

import OpenCL

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
