//
//  Grid.swift
//
//
//  Created by Philip Turner on 5/16/24.
//

/// A series of fine uniform grids spanning one cubic Bohr.
public struct Voxel {
  /// A list of levels, going from coarse to the finest level available.
  public var levels: [Level] = []
  
  public init() {
    
  }
}

/// A coarse, uniform grid encapsulating a region of the domain.
public struct Grid {
  /// A list of levels, going from fine (1x1x1 Bohr) to the coarsest level
  /// available.
  public var levels: [Level] = []
  
  /// The remaining levels of the grid, which are allocated sparsely.
  ///
  /// The dimensions must match the 1x1x1 level.
  public var voxels: [Voxel?] = []
  
  public init() {
    
  }
}
