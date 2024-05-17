//
//  Voxel.swift
//
//
//  Created by Philip Turner on 5/16/24.
//

struct VoxelDescriptor {
  // The number of mipmap levels.
  var depth: Int?
}

// A data structure encapsulating one cubic Bohr.
struct Voxel {
  // A list of levels, going from coarse to the finest level available.
  var levels: [Level] = []
  
  init(descriptor: VoxelDescriptor) {
    guard let depth = descriptor.depth else {
      fatalError("Descriptor was invalid.")
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
