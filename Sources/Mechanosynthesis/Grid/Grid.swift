//
//  Grid.swift
//
//
//  Created by Philip Turner on 5/16/24.
//

public struct GridDescriptor {
  /// The start of the smallest bounding box that encloses the data.
  ///
  /// Units: 1 cubic Bohr
  public var offset: SIMD3<Int32>?
  
  /// The size of the smallest bounding box that encloses the data.
  ///
  /// Units: 1 cubic Bohr
  public var dimensions: SIMD3<UInt32>?
  
  public init() {
    
  }
}

/// A series of fine uniform grids spanning one cubic Bohr.
public struct Voxel {
  /// A list of levels, going from coarse to the finest level available.
  public var levels: [Level] = []
  
  public init() {
    
  }
}

/// A coarse, uniform grid encapsulating a region of the domain.
///
/// The offset and dimensions are specified in integer multiples of the Bohr
/// radius.
public struct Grid {
  /// The start of the smallest bounding box that encloses the data.
  public var offset: SIMD3<Int32>
  
  /// The size of the smallest bounding box that encloses the data.
  public var dimensions: SIMD3<UInt32>
  
  /// Encapsulates the data for valid 1x1x1 Bohr cells.
  public var highestLevel: Level
  
  /// The remaining levels of the grid, which are allocated sparsely.
  public var voxels: [Voxel?]
  
  public init(descriptor: GridDescriptor) {
    guard let offset = descriptor.offset,
          let dimensions = descriptor.dimensions else {
      fatalError("Descriptor was incomplete.")
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
    let voxelCount = Int(dimensions[0] * dimensions[1] * dimensions[2])
    voxels = Array(repeating: nil, count: voxelCount)
  }
}
