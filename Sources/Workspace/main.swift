//
//  main.swift
//
//
//  Created by Philip Turner on 5/25/24.
//

import Mechanosynthesis
import Numerics

// Tasks:
// - Create code that can fill in the wavefunction amplitude at each voxel.
// - Test out OrbitalMesh on the origin-centered N atom.
// - Create code that can shift the nuclei into a mostly centered coordinate
//   space, creating a good octree. Then, fuse the shift into 'grid.offset'.
// - Test out OrbitalMesh on the N2 molecule.
// - Create code that can fuse the meshes from a large number of electrons into
//   a global simulation domain.
//   - There would be a "master" mesh whose amplitude is zero everywhere.

// Create an ansatz for the electronic structure.
var ansatzDesc = AnsatzDescriptor()
ansatzDesc.atomicNumber = 1
ansatzDesc.netSpin = 0.5
ansatzDesc.position = .zero
ansatzDesc.sizeExponent = 4
let ansatz = Ansatz(descriptor: ansatzDesc)

let orbital = ansatz.orbitals[0]
var mesh = OrbitalMesh(orbital: orbital)

// Fill in the voxels.
for nodeID in orbital.octree.nodes.indices {
  let node = orbital.octree.nodes[nodeID]
  guard node.spacing == 1 else {
    // Only consider nodes where the spacing is 1 Bohr.
    continue
  }
  
  // Locate this chunk within the grid.
  let voxelOffsetFloat = node.center - SIMD3<Float>(mesh.grid.offset)
  let voxelOffset = SIMD3<Int>(voxelOffsetFloat.rounded(.down))
  
  var voxelLinearIndex: Int = .zero
  do {
    let dimensions = mesh.grid.dimensions
    voxelLinearIndex += voxelOffset[0]
    voxelLinearIndex += voxelOffset[1] * dimensions[0]
    voxelLinearIndex += voxelOffset[2] * dimensions[0] * dimensions[1]
  }
  
  // Create an empty voxel.
  var voxel = Voxel()
  
  // Search through the children recursively.
  func traverseOctreeNode(nodeID: UInt32, levelID: UInt32) {
    // Ensure the current level is allocated.
    if voxel.levels.count == levelID {
      var levelDesc = LevelDescriptor()
      let gridSize = 1 << levelID
      levelDesc.dimensions = SIMD3(repeating: gridSize)
      let level = Level(descriptor: levelDesc)
      voxel.levels.append(level)
      
      print("Allocated a new level: \(levelID)")
    }
    
    // Retrieve the node.
    let node = mesh.orbital.octree.nodes[Int(nodeID)]
    print("\(levelID) - \(node.center)")
    
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
    
    // Iterate over the children.
    for branchID in 0..<8 {
      let childOffset = node.branchesMask[branchID]
      guard childOffset < 255 else {
        continue
      }
      let childOffset32 = UInt32(truncatingIfNeeded: childOffset)
      let childNodeID = node.branchesIndex + childOffset32
      print(nodeID, childNodeID)
      
      traverseOctreeNode(nodeID: childNodeID, levelID: levelID + 1)
    }
  }
  
  print()
  traverseOctreeNode(nodeID: UInt32(nodeID), levelID: 0)
}
