//
//  OrbitalMesh.swift
//
//
//  Created by Philip Turner on 5/25/24.
//

import Mechanosynthesis

// This does not belong in the 'Mechanosynthesis' library. It should be copied
// from source whenever it is used.
struct OrbitalMesh {
  var orbital: HydrogenicOrbital
  var grid: Grid
  
  init(orbital: HydrogenicOrbital) {
    self.orbital = orbital
    self.grid = Self.createGrid(orbital: orbital)
    
    createHighestLevel()
    createVoxels()
  }
  
  // Creates an empty grid with the required bounds.
  static func createGrid(orbital: HydrogenicOrbital) -> Grid {
    // Allocate variables to accumulate the mesh bounds.
    var minimumBound: SIMD3<Int32> = .init(repeating: .max)
    var maximumBound: SIMD3<Int32> = .init(repeating: -.max)

    // Iterate over the octree nodes.
    let orbital = ansatz.orbitals[0]
    for node in orbital.octree.nodes {
      guard node.spacing == 2 else {
        // Only consider nodes where the spacing is 2 Bohr.
        continue
      }
      
      // Find the bounds of this node.
      let lowerCorner = SIMD3<Int32>(node.center - node.spacing / 2)
      let upperCorner = SIMD3<Int32>(node.center + node.spacing / 2)
      
      // Merge with the bounds of the entire mesh.
      minimumBound
        .replace(with: lowerCorner, where: lowerCorner .< minimumBound)
      maximumBound
        .replace(with: upperCorner, where: upperCorner .> maximumBound)
    }
    
    // Create the grid.
    var gridDesc = GridDescriptor()
    gridDesc.offset = minimumBound
    gridDesc.dimensions = SIMD3<UInt32>(
      truncatingIfNeeded: maximumBound &- minimumBound)
    return Grid(descriptor: gridDesc)
  }
  
  // Fill in the highest level of the grid.
  mutating func createHighestLevel() {
    for node in orbital.octree.nodes {
      guard node.spacing == 2 else {
        // Only consider nodes where the spacing is 2 Bohr.
        continue
      }
      
      // Locate this chunk within the grid.
      let lowerCorner = SIMD3<Int32>(node.center - node.spacing / 2)
      let voxelIndexOffset = SIMD3<UInt32>(
        truncatingIfNeeded: lowerCorner &- grid.offset)
      let chunkIndex = voxelIndexOffset / 2
      
      var chunkLinearIndex: UInt32 = .zero
      do {
        let dimensions = grid.dimensions / 2
        chunkLinearIndex += chunkIndex[0]
        chunkLinearIndex += chunkIndex[1] * dimensions[0]
        chunkLinearIndex += chunkIndex[2] * dimensions[0] * dimensions[1]
      }
      
      // Calculate the wavefunction amplitude for each child.
      var x = SIMD8<Float>(0, 1, 0, 1, 0, 1, 0, 1) * 0.5 - 0.25
      var y = SIMD8<Float>(0, 0, 1, 1, 0, 0, 1, 1) * 0.5 - 0.25
      var z = SIMD8<Float>(0, 0, 0, 0, 1, 1, 1, 1) * 0.5 - 0.25
      x = x * node.spacing + node.center.x
      y = y * node.spacing + node.center.y
      z = z * node.spacing + node.center.z
      var amplitude = orbital.basisFunction.amplitude(x: x, y: y, z: z)
      
      // Mark unoccupied cells with NAN.
      let mask32 = SIMD8<UInt32>(truncatingIfNeeded: node.branchesMask)
      amplitude.replace(with: .nan, where: mask32 .!= 255)
      
      // Write data for every cell that terminates at 1x1x1.
      grid.highestLevel.data[Int(chunkLinearIndex)] = amplitude
    }
  }
  
  // Fill in the multigrid for each voxel.
  mutating func createVoxels() {
    for nodeID in orbital.octree.nodes.indices {
      let node = orbital.octree.nodes[nodeID]
      guard node.spacing == 1 else {
        // Only consider nodes where the spacing is 1 Bohr.
        continue
      }
      
      // Create an empty voxel.
      var voxel = Voxel()
      let voxelMinimum = node.center.rounded(.down)
      
      // Search through the children recursively.
      func traverseOctreeNode(nodeID: UInt32, levelID: UInt32) {
        // Ensure the current level is allocated.
        if voxel.levels.count == levelID {
          var levelDesc = LevelDescriptor()
          let gridSize = UInt32(1) << levelID
          levelDesc.dimensions = SIMD3(repeating: gridSize)
          let level = Level(descriptor: levelDesc)
          voxel.levels.append(level)
        }
        
        // Retrieve the node.
        let node = orbital.octree.nodes[Int(nodeID)]
        
        // Iterate over the children.
        for branchID in 0..<8 {
          let childOffset = node.branchesMask[branchID]
          guard childOffset < 255 else {
            continue
          }
          let childOffset32 = UInt32(truncatingIfNeeded: childOffset)
          let childNodeID = node.branchesIndex + childOffset32
          traverseOctreeNode(nodeID: childNodeID, levelID: levelID + 1)
        }
        
        // Calculate the wavefunction amplitude for each child.
        var x = SIMD8<Float>(0, 1, 0, 1, 0, 1, 0, 1) * 0.5 - 0.25
        var y = SIMD8<Float>(0, 0, 1, 1, 0, 0, 1, 1) * 0.5 - 0.25
        var z = SIMD8<Float>(0, 0, 0, 0, 1, 1, 1, 1) * 0.5 - 0.25
        x = x * node.spacing + node.center.x
        y = y * node.spacing + node.center.y
        z = z * node.spacing + node.center.z
        var amplitude = orbital.basisFunction.amplitude(x: x, y: y, z: z)
        
        // Mark unoccupied cells with NAN.
        let mask32 = SIMD8<UInt32>(truncatingIfNeeded: node.branchesMask)
        amplitude.replace(with: .nan, where: mask32 .!= 255)
        
        // Locate this node within the level.
        var nodeLinearIndex: UInt32 = .zero
        do {
          var nodeOffsetFloat = node.center - voxelMinimum
          let gridSize = UInt32(1) << levelID
          nodeOffsetFloat *= Float(gridSize)
          
          let nodeOffset = SIMD3<UInt32>(nodeOffsetFloat.rounded(.down))
          let dimensions = SIMD3<UInt32>(repeating: gridSize)
          nodeLinearIndex += nodeOffset[0]
          nodeLinearIndex += nodeOffset[1] * dimensions[0]
          nodeLinearIndex += nodeOffset[2] * dimensions[0] * dimensions[1]
        }
        
        // Write data for every cell that terminates at this level.
        // - Doing this in a single statement, to avoid memory copies from CoW.
        voxel.levels[Int(levelID)]
          .data[Int(nodeLinearIndex)] = amplitude
      }
      traverseOctreeNode(nodeID: UInt32(nodeID), levelID: 0)
      
      // Locate this chunk within the grid.
      var voxelLinearIndex: UInt32 = .zero
      do {
        let voxelOffsetFloat = node.center - SIMD3<Float>(grid.offset)
        let voxelOffset = SIMD3<UInt32>(voxelOffsetFloat.rounded(.down))
        let dimensions = grid.dimensions
        voxelLinearIndex += voxelOffset[0]
        voxelLinearIndex += voxelOffset[1] * dimensions[0]
        voxelLinearIndex += voxelOffset[2] * dimensions[0] * dimensions[1]
      }
      grid.voxels[Int(voxelLinearIndex)] = voxel
    }
  }
}
