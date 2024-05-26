//
//  Level.swift
//
//
//  Created by Philip Turner on 5/16/24.
//

public struct LevelDescriptor {
  /// The start of the smallest bounding box that encloses the data.
  ///
  /// Units: cell spacing (twice as fine as chunk spacing)
  public var offset: SIMD3<Int32>?
  
  /// The size of the smallest bounding box that encloses the data.
  ///
  /// Units: cell spacing (twice as fine as chunk spacing)
  public var dimensions: SIMD3<UInt32>?
  
  public init() {
    
  }
}

/// A uniform grid encapsulating one mipmap level of a voxel.
public struct Level {
  /// The start of the smallest bounding box that encloses the data.
  public var offset: SIMD3<Int32>
  
  /// The size of the smallest bounding box that encloses the data.
  public var dimensions: SIMD3<UInt32>
  
  /// The chunks in the level.
  ///
  /// Reorders data at the 2x2x2 granularity, to improve memory locality and
  /// decrease the overhead of dispatching compute work. The cells within
  /// each 2x2x2 chunk are stored in Morton order.
  ///
  /// Unoccupied cells have `NAN` for the data value.
  public var data: [SIMD8<Float>]
  
  public init(descriptor: LevelDescriptor) {
    guard let offset = descriptor.offset,
          let dimensions = descriptor.dimensions else {
      fatalError("Descriptor was incomplete.")
    }
    guard all(dimensions .> 0) else {
      fatalError("Chunk count must be nonzero.")
    }
    guard all(offset % 2 .== 0),
          all(dimensions % 2 .== 0) else {
      fatalError("Level must be aligned to a multiple of 2x2x2 cells.")
    }
    self.offset = offset
    self.dimensions = dimensions
    
    // Allocate an array of chunks.
    let chunkDimensions = dimensions / 2
    let chunkCount = Int(
      chunkDimensions[0] * chunkDimensions[1] * chunkDimensions[2])
    data = Array(
      repeating: SIMD8(repeating: .nan), count: chunkCount)
  }
}
