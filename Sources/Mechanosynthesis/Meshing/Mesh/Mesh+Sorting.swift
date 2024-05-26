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
        let mappedLocation = SIMD2<UInt32>(repeating: .max)
        nodesToVoxelsMap.append(mappedLocation)
        continue
      }
      
      // Locate this node within the grid of voxels.
      let coordinate = node.center / Float(spacing)
      let coordinateDelta = coordinate - SIMD3<Float>(coarseVoxels.offset)
      guard all(coordinateDelta .> 0),
            all(coordinateDelta .< SIMD3(coarseVoxels.dimensions)) else {
        fatalError("Voxel was out of bounds: \(coordinate) \(coarseVoxels.offset).")
      }
      var voxelID: UInt32 = .zero
      do {
        let deltaInt = SIMD3<UInt32>(coordinateDelta.rounded(.down))
        let dimensions = coarseVoxels.dimensions
        voxelID += deltaInt[0]
        voxelID += deltaInt[1] * dimensions[0]
        voxelID += deltaInt[2] * dimensions[0] * dimensions[1]
      }
      
      // Accumulate the number of nodes mapped to this voxel.
      let slotID = voxelAccumulators[Int(voxelID)]
      voxelAccumulators[Int(voxelID)] = slotID + 1
      let mappedLocation = SIMD2(voxelID, slotID)
      nodesToVoxelsMap.append(mappedLocation)
    }
    print(voxelAccumulators)
    print(nodesToVoxelsMap.count)
    print(octree.nodes.count)
    
    fatalError("Not implemented.")
  }
}
