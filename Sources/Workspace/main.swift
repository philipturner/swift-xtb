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
for node in orbital.octree.nodes {
  guard node.spacing <= 1 else {
    // Only consider nodes where the spacing is 1 Bohr or less.
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
  
  // If the voxel does not already exist, create it.
  if mesh.grid.voxels[voxelLinearIndex] == nil {
    
    let voxel = Voxel()
    mesh.grid.voxels[voxelLinearIndex] = voxel
  }
  var voxel = mesh.grid.voxels[voxelLinearIndex]
  
  // TODO: Respect copy-on-write semantics when the voxels get really deep,
  // and touching them the wrong way incurs a large memory copy.
  
  // Allocate memory for the wavefunction amplitude data.
}
