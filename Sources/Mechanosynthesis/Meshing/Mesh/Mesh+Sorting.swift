//
//  Mesh+Sorting.swift
//
//
//  Created by Philip Turner on 5/26/24.
//

extension Mesh {
  // Maps octree nodes to coarse voxels.
  //
  // Unsure how to proceed with fine voxels for now.
  func mapNodesToVoxels(octree: Octree) -> [SIMD2<UInt32>] {
    // Create an accumulator for each voxel, just for this octree.
    var voxelAccumulators = [UInt32](
      repeating: .zero, count: coarseVoxels.cells.count)
    
    // Iterate over the nodes.
    var nodesToVoxelsMap: [SIMD2<UInt32>] = []
    for node in octree.nodes {
      guard node.spacing <= Float(spacing) else {
        continue
      }
      
      // Locate this node within the grid of voxels.
      let coordinate = node.center / Float(spacing)
      let coordinateDelta = coordinate - SIMD3<Float>(coarseVoxels.offset)
      guard all(coordinateDelta .> 0),
            all(coordinateDelta .< SIMD3(coarseVoxels.dimensions)) else {
        fatalError("Voxel was out of bounds: \(coordinate) \(coarseVoxels.offset).")
      }
      print("success")
    }
    
    fatalError("Not implemented.")
  }
}
