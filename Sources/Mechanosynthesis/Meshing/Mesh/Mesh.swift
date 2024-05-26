//
//  Mesh.swift
//
//
//  Created by Philip Turner on 5/26/24.
//

public struct MeshDescriptor {
  // Required. The set of octrees to generate the topology from.
  public var octrees: [Octree] = []
  
  // Required. The power-2 size of coarse voxels, in Bohr.
  public var sizeExponent: Int?
  
  public init() {
    
  }
}

public struct Mesh {
  // A uniform grid of coarse voxels.
  public var coarseVoxels: Grid<CoarseVoxel>
  
  // The separation between coarse voxels, in Bohr.
  public var spacing: Int
  
  public init(descriptor: MeshDescriptor) {
    fatalError("Not implemented.")
  }
}

extension Mesh {
  public mutating func append(contentsOf other: Mesh) {
    fatalError("Not implemented.")
  }
}
