import XCTest
import Mechanosynthesis
import Numerics

// Ideas for test systems:
// - F2 molecule with an F* *F ansatz
//   - Span the longest diagonal of a fine voxel, to maximize the number of
//     performance edge cases triggered.
// - C2H2 molecule with an H* *C* *C* *H ansatz
//   - Omit orbitals with m > -l, ensure results are the same. Establish this
//     as good practice, but it's something the user should do. The library
//     should not be responsible for de-duplicating octrees.
//
// Ideas for test coverage:
// - Depict a cross-section of each resolution level's span in comments.
//   - Create a function to automate the rendering.
//   - Examine how well different coarse voxel sizes conform to the
//     non-cuboidal structure of valence orbitals.
// - Quantifying the number of cells at each resolution level, ensuring the
//   same number is always produced.
// - Examining the octree generation vs. mesh generation time.
