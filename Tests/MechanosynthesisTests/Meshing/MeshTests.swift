import XCTest
import Mechanosynthesis
import Numerics

// Ideas for test systems:
// - F2 molecule with an F* *F ansatz
//   - Span the longest diagonal of a fine voxel, to maximize the number of
//     performance edge cases triggered.
// - C2H2 molecule with an H* *C* *C* *H ansatz.
//   - Study how the rotational asymmetry of the p-orbitals impacts the
//     generated mesh. Superimpose three rotated versions of the octree for
//     spatially polarized orbitals.
//   - Provide an API to control the orientation of cubic harmonics in Ansatz
//     (principal axes; Float32, unlike MM4RigidBody).
//
// Ideas for test coverage:
// - Depict a cross-section of each resolution level's span in comments.
//   - Create a function to automate the rendering.
//   - Examine how well different coarse voxel sizes conform to the
//     non-cuboidal structure of valence orbitals.
// - Quantifying the number of cells at each resolution level, ensuring the
//   same number is always produced.
// - Examining the octree generation vs. mesh generation time.
