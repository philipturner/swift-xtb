//
//  OrbitalMesh.swift
//
//
//  Created by Philip Turner on 5/25/24.
//

import Mechanosynthesis

struct OrbitalMesh {
  var orbital: HydrogenicOrbital
  var grid: Grid
  
  init(orbital: HydrogenicOrbital) {
    self.orbital = orbital
    self.grid = Self.createGrid(orbital: orbital)
    
    fillHighestLevel()
  }
  
  // Creates an empty grid with the required bounds.
  static func createGrid(orbital: HydrogenicOrbital) -> Grid {
    // Allocate variables to accumulate the mesh bounds.
    var minimumBound: SIMD3<Int> = .init(repeating: .max)
    var maximumBound: SIMD3<Int> = .init(repeating: -.max)

    // Iterate over the octree nodes.
    let orbital = ansatz.orbitals[0]
    for node in orbital.octree.nodes {
      guard node.spacing == 2 else {
        // Only consider nodes where the spacing is 2 Bohr.
        continue
      }
      
      // Find the bounds of this node.
      let lowerCorner = SIMD3<Int>(node.center - node.spacing / 2)
      let upperCorner = SIMD3<Int>(node.center + node.spacing / 2)
      
      // Merge with the bounds of the entire mesh.
      minimumBound.replace(with: lowerCorner, where: lowerCorner .< minimumBound)
      maximumBound.replace(with: upperCorner, where: upperCorner .> maximumBound)
    }
    
    // Create the grid.
    var gridDesc = GridDescriptor()
    gridDesc.offset = minimumBound
    gridDesc.dimensions = maximumBound &- minimumBound
    return Grid(descriptor: gridDesc)
  }
  
  // Fill in the highest level of the grid.
  mutating func fillHighestLevel() {
    for node in orbital.octree.nodes {
      guard node.spacing == 2 else {
        // Only consider nodes where the spacing is 2 Bohr.
        continue
      }
      
      // Locate this chunk within the grid.
      let lowerCorner = SIMD3<Int>(node.center - node.spacing / 2)
      let voxelIndexOffset = lowerCorner &- grid.offset
      let chunkIndex = voxelIndexOffset / 2
      
      var chunkLinearIndex: Int = .zero
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
      grid.highestLevel.data[chunkLinearIndex] = amplitude
    }
  }
}
