//
//  main.swift
//
//
//  Created by Philip Turner on 5/25/24.
//

import Mechanosynthesis
import Numerics

// Tasks:
// - Minimize the bottleneck from deep levels allocating a lot of memory in the
//   sparse voxels.
// - Create code that can shift the nuclei into a mostly centered coordinate
//   space, creating a good octree. Then, fuse the shift into 'grid.offset'.
// - Test out AnsatzMesh on the N2 molecule.
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

// Report the number of cells in the entire mesh. Has it decreased from 300k
// for the 1s orbital of nitrogen?
