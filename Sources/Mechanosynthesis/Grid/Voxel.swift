//
//  Voxel.swift
//
//
//  Created by Philip Turner on 5/16/24.
//

public struct VoxelDescriptor {
  /// The number of mipmap levels.
  public var depth: Int?
  
  public init() {
    
  }
}

/// A series of fine uniform grids spanning one cubic Bohr.
public struct Voxel {
  /// A list of levels, going from coarse to the finest level available.
  public var levels: [Level] = []
  
  public init(descriptor: VoxelDescriptor) {
    guard let depth = descriptor.depth else {
      fatalError("Descriptor was incomplete.")
    }
    guard depth > 0 else {
      fatalError("Voxel must contain at least one level.")
    }
    
    // Allocate an array of levels.
    for levelID in 0..<depth {
      let dimensions = SIMD3<Int>(repeating: 1 << levelID)
      
      var levelDesc = LevelDescriptor()
      levelDesc.dimensions = dimensions
      let level = Level(descriptor: levelDesc)
      levels.append(level)
    }
  }
}
