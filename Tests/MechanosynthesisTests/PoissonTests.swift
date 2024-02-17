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
