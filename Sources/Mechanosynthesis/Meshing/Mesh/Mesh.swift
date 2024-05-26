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
    guard let sizeExponent = descriptor.sizeExponent else {
      fatalError("Descriptor was incomplete.")
    }
    guard sizeExponent > 0 else {
      fatalError("Coarse voxel spacing must be at least 2 Bohr.")
    }
    spacing = 1 << sizeExponent
    
    // Find the bounding box of the coarse grid.
    print()
    var minimum: SIMD3<Int32> = .init(repeating: .max)
    var maximum: SIMD3<Int32> = .init(repeating: -.max)
    for octree in descriptor.octrees {
      // The highest octree level must span 2 * coarse voxel spacing,
      // so that either half of the octree is 1 * coarse voxel spacing.
      let expectedSpacing = Float(2 * spacing)
      guard octree.nodes[0].spacing >= expectedSpacing else {
        fatalError("Octree must be large enough to fit the coarse voxels.")
      }
      
      for node in octree.nodes {
        // There are some topologies where no nodes will meet this criterion.
        // For example, if the octree terminates with only a single node.
        //
        // Remember, the coarse voxel spacing is twice the resolution of the
        // actual data that is stored.
        guard node.spacing == Float(spacing) else {
          continue
        }
        
        // Find the bounding box of the node.
        let nodeMinimum = SIMD3<Int32>(node.center - node.spacing / 2)
        let nodeMaximum = SIMD3<Int32>(node.center + node.spacing / 2)
        minimum.replace(with: nodeMinimum, where: nodeMinimum .< minimum)
        maximum.replace(with: nodeMaximum, where: nodeMaximum .> maximum)
        print("node:")
        print("\(node.center - node.spacing / 2) -> \(nodeMinimum)")
        print("\(node.center + node.spacing / 2) -> \(nodeMaximum)")
      }
    }
    guard all(minimum .< maximum) else {
      // This happens when there are no octrees.
      fatalError("Mesh bounds could not be established.")
    }
    print()
    print("mesh:")
    print("minimum: \(minimum)")
    print("maximum: \(maximum)")
    
    fatalError("Not implemented.")
  }
}

extension Mesh {
  public mutating func append(contentsOf other: Mesh) {
    fatalError("Not implemented.")
  }
}
