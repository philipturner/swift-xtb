import XCTest
import Mechanosynthesis
import Numerics

// Prove that a multigrid solver with rough mipmap sampling works. This is
// different from the Real-Space DFT paper, which used a complicated
// interpolation scheme.
//
// The combined electron and ion Hartree potential can be used to compute atomic
// forces. Interpret the potential as a voltage, and the derivative of voltage
// is electric field. The electric field exerts a force on the nucleus. The
// cost is O(1) * O(number of nuclei).
//
// 64x64x64 with the current multigrid algorithm idea, by default.
// 65x65x65 for debugging with the Real-Space DFT paper algorithm, if
// necessary.
final class PoissonTests {
  // Truncate at the quadrupole moment because it's the last one frequently
  // given in analytical form. Set the origin to the "center of charge" to
  // have defined behavior when the quadrupole depends on the origin.
  // https://phys.libretexts.org/Bookshelves/Mathematical_Physics_and_Pedagogy/Mathematical_Methods/The_Multipole_Expansion
  func testMultipoleExpansion() throws {
    // Create a list of point charges, compute the multipole expansion, and
    // compare to results from direct integration. Show how the accuracy
    // improves from cutoff -> monopole -> dipole -> quadrupole.
    
  }
}
