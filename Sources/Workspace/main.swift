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
for voxel in mesh.grid.voxels where voxel != nil {
  print(voxel!.levels.count)
}
