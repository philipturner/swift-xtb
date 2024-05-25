//
//  Level.swift
//
//
//  Created by Philip Turner on 5/16/24.
//

public struct LevelDescriptor {
  // The number of chunks along each axis.
  public var dimensions: SIMD3<Int>?
  
  public init() {
    
  }
}

/// A uniform grid encapsulating one mipmap level of a voxel.
public struct Level {
  /// The chunks in the level.
  ///
  /// Reorders data at the 2x2x2 granularity, to improve memory locality and
  /// decrease the overhead of dispatching compute work. The cells within
  /// each 2x2x2 chunk are stored in Morton order.
  public var data: [SIMD8<Float>]
  
  /// Whether each cell is part of the "valid" region.
  public var mask: [UInt8]
  
  public init(descriptor: LevelDescriptor) {
    guard let dimensions = descriptor.dimensions else {
      fatalError("Descriptor was incomplete.")
    }
    guard all(dimensions .> 0) else {
      fatalError("Chunk count must be nonzero.")
    }
    
    // Allocate an array of chunks.
    let chunkCount = dimensions[0] * dimensions[1] * dimensions[2]
    data = Array(repeating: SIMD8(repeating: .nan), count: chunkCount)
    mask = Array(repeating: 0b0000_0000, count: chunkCount)
  }
}
