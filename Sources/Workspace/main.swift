//
//  main.swift
//
//
//  Created by Philip Turner on 5/25/24.
//

import Mechanosynthesis
import Numerics

// Tasks:
// - Draft the code for the new "HydrogenicOrbitalMesh" API.
//   - Avoid copying/saving the wavefunction amplitudes this time. The most
//     important information is just zeroes and NANs, stored as Float32.
// - Test out HydrogenicOrbitalMesh on the N atom.
//   - Report the number of cells in the entire mesh. Has it decreased from
//     300k for the 1s orbital?
// - Create code that can shift the nuclei into a mostly centered coordinate
//   space, creating a good octree. Then, fuse the shift into 'grid.offset'.
// - Test out AnsatzMesh on the N2 molecule.
// - Create code that can fuse the meshes from a large number of electrons into
//   a global simulation domain.
//   - There would be a "master" mesh whose amplitude is zero everywhere.
