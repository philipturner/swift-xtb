import XCTest
import Mechanosynthesis
import Numerics

// Ideas for test systems:
// - N2 molecule with a (N-)** **(N+) ansatz.
// - H2O molecule with a H* *O* *H ansatz.
// - C2H2 molecule with a H* *(C-) (C+)* *H ansatz.
//
// Ideas for test coverage:
// - Depict a cross-section of each resolution level's span in comments.
//   - Create a function to automate the rendering.
//   - Examine how well different coarse voxel sizes conform to the
//     non-cuboidal structure of valence orbitals.
// - Quantifying the number of cells at each resolution level, ensuring the
//   same number is always produced.
// - Examining the octree generation vs. mesh generation time.
