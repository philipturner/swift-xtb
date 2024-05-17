//
//  Grid.swift
//
//
//  Created by Philip Turner on 5/16/24.
//

struct GridDescriptor {
  // The start of the smallest bounding box that encloses the data.
  //
  // Units: 1 cubic Bohr
  var offset: SIMD3<Int>?
  
  // The size of the smallest bounding box that encloses the data.
  //
  // Units: 1 cubic Bohr
  var dimensions: SIMD3<Int>?
}

// A uniform grid encapsulating a compact region of the domain.
//
// The offset and dimensions are specified in integer multiples of the Bohr
// radius.
struct Grid {
  // The start of the smallest bounding box that encloses the data.
  var offset: SIMD3<Int>
  
  // The size of the smallest bounding box that encloses the data.
  var dimensions: SIMD3<Int>
  
  // Encapsulates the data for valid 1x1x1 Bohr cells.
  var highestLevel: Level
  
  // The remaining levels of the grid, which are allocated sparsely.
  var voxels: [Voxel?]
  
  init(descriptor: GridDescriptor) {
    guard let offset = descriptor.offset,
          let dimensions = descriptor.dimensions else {
      fatalError("Descriptor was invalid.")
    }
    guard all(dimensions .> 0) else {
      fatalError("Voxel count must be nonzero.")
    }
    guard all(offset % 2 .== 0),
          all(dimensions % 2 .== 0) else {
      fatalError("Grid must be aligned to a multiple of 2x2x2 Bohr.")
    }
    self.offset = offset
    self.dimensions = dimensions
    
    // Allocate a level for the coarsest data.
    var levelDesc = LevelDescriptor()
    levelDesc.dimensions = dimensions / 2
    highestLevel = Level(descriptor: levelDesc)
    
    // Allocate an array of voxels.
    let voxelCount = dimensions[0] * dimensions[1] * dimensions[2]
    voxels = Array(repeating: nil, count: voxelCount)
  }
}
