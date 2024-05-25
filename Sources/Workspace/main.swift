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
  print(minimumBound, maximumBound)
}

// Initialize an empty grid.
var gridDesc = GridDescriptor()
gridDesc.offset = minimumBound
gridDesc.dimensions = maximumBound &- minimumBound
let grid = Grid(descriptor: gridDesc)
print(grid.voxels.count)
