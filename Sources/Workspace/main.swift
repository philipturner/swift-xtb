//
//  main.swift
//
//
//  Created by Philip Turner on 5/25/24.
//

import Mechanosynthesis
import Numerics

// Tasks:
// - Learn how to set up a mesh for the hydrogen atom.
// - Wrap the grid generation code into a function, store in 'Workspace'.

// Create an ansatz for the electronic structure.
var ansatzDesc = AnsatzDescriptor()
ansatzDesc.atomicNumber = 1
ansatzDesc.netSpin = 0.5
ansatzDesc.position = .zero
ansatzDesc.sizeExponent = 4
let ansatz = Ansatz(descriptor: ansatzDesc)

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

// Initialize an empty grid.
var gridDesc = GridDescriptor()
gridDesc.offset = minimumBound
gridDesc.dimensions = maximumBound &- minimumBound
var grid = Grid(descriptor: gridDesc)

// Fill in the highest level of the grid.
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
  let chunkIndexBounds = grid.dimensions / 2
  chunkLinearIndex += chunkIndex[0]
  chunkLinearIndex += chunkIndex[1] * chunkIndexBounds[0]
  chunkLinearIndex += chunkIndex[2] * chunkIndexBounds[0] * chunkIndexBounds[1]
  
  // Calculate the wavefunction amplitude for each child.
  var x = SIMD8<Float>(0, 1, 0, 1, 0, 1, 0, 1) * 0.5 - 0.25
  var y = SIMD8<Float>(0, 0, 1, 1, 0, 0, 1, 1) * 0.5 - 0.25
  var z = SIMD8<Float>(0, 0, 0, 0, 1, 1, 1, 1) * 0.5 - 0.25
  x = node.spacing * x + node.center[0]
  y = node.spacing * y + node.center[1]
  z = node.spacing * z + node.center[2]
  var amplitude = orbital.basisFunction.amplitude(x: x, y: y, z: z)
  
  // Mark unoccupied cells with NAN.
  let mask32 = SIMD8<UInt32>(truncatingIfNeeded: node.branchesMask)
  amplitude.replace(with: .nan, where: mask32 .!= 255)
  
  // Write data for every cell that terminates at 1x1x1.
  grid.highestLevel.data[chunkLinearIndex] = amplitude
}

// Fill in the voxels.
