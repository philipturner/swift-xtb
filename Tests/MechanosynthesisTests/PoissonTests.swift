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
final class PoissonTests: XCTestCase {
  // Include at least the quadrupole expansion in the final code. Set the
  // origin to the "center of charge" to have defined behavior when the
  // quadrupole depends on the origin.
  // https://phys.libretexts.org/Bookshelves/Mathematical_Physics_and_Pedagogy/Mathematical_Methods/The_Multipole_Expansion
  //
  // Also reproducing the octupole expansion. My instinct tells me that this
  // could be a better performance tradeoff than quadrupoles. However, the
  // decision of where to truncate the expansion should follow a rigorous
  // investigation.
  // https://physics.stackexchange.com/questions/269753/electric-octupole-moment-in-cartesian-coordinates
  
  // Hexadecapole expansion on this page, formula 19:
  // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7017380/
  // Will not be implementing this due to the high complexity of the source
  // code.
  func testMultipoleExpansion() throws {
    // Create a list of point charges, compute the multipole expansion, and
    // compare to results from direct integration. Show how the accuracy
    // improves from cutoff -> monopole -> dipole -> quadrupole -> direct.
    
    // Program settings.
    let ionCharge: Int = 0
    let samplePosition: SIMD3<Float> = [0, 0, 2]
    
    // Create the list of charges.
    var pointCharges: [SIMD4<Float>] = []
    pointCharges.append(SIMD4(-0.1, 0.2, 0.2, 7))
    for _ in 0..<2 {
      pointCharges.append(SIMD4(0.1, 0.2, 0, -1))
      pointCharges.append(SIMD4(-0.3, 0, -0.5, -1))
      pointCharges.append(SIMD4(-0.01, 0.4, 0.03, -1))
    }
    for _ in 0..<(1 - ionCharge) {
      pointCharges.append(SIMD4(0.3, 0.4, 0.5, -1))
    }
    
    // Find the origin of the multipole expansion.
    var origin: SIMD3<Float> = .zero
    var totalMass: Float = .zero
    for charge in pointCharges {
      let position = unsafeBitCast(charge, to: SIMD3<Float>.self)
      let mass = charge.w.magnitude
      origin += position * mass
      totalMass += mass
    }
    origin /= totalMass
    
    // Pre-compute some quantities.
    let r = samplePosition - origin
    let rNorm = (r * r).sum().squareRoot()
    let rHat = r / rNorm
    
    // Record each estimate to the electric potential, starting with 'cutoff'.
    var estimates: [Float] = []
    var potential: Float = .zero
    estimates.append(potential)
    
    // Compute the monopole expansion.
    do {
      var totalCharge: Float = .zero
      for charge in pointCharges {
        totalCharge += charge.w
      }
      potential += totalCharge / rNorm
    }
    estimates.append(potential)
    
    // Compute the dipole expansion.
    do {
      var dipoleMoment: SIMD3<Float> = .zero
      for charge in pointCharges {
        let position = unsafeBitCast(charge, to: SIMD3<Float>.self)
        let rPrime = position - origin
        dipoleMoment += charge.w * rPrime
      }
      let dotProduct = (rHat * dipoleMoment).sum()
      potential += dotProduct / (rNorm * rNorm)
    }
    estimates.append(potential)
    
    // Compute the quadrupole expansion.
    do {
      var quadrupoleMoment: (SIMD3<Float>, SIMD3<Float>, SIMD3<Float>)
      quadrupoleMoment = (.zero, .zero, .zero)
      
      for charge in pointCharges {
        let position = unsafeBitCast(charge, to: SIMD3<Float>.self)
        let rPrime = position - origin
        let rPrimeNorm = (rPrime * rPrime).sum().squareRoot()
        
        // Activate the trace part after debugging what you have here.
        var tensor = (
          3 * rPrime.x * rPrime,
          3 * rPrime.y * rPrime,
          3 * rPrime.z * rPrime)
        tensor.0[0] -= 1 * rPrimeNorm * rPrimeNorm
        tensor.1[1] -= 1 * rPrimeNorm * rPrimeNorm
        tensor.2[2] -= 1 * rPrimeNorm * rPrimeNorm
        
        quadrupoleMoment.0 += charge.w * tensor.0
        quadrupoleMoment.1 += charge.w * tensor.1
        quadrupoleMoment.2 += charge.w * tensor.2
      }
      
      let QrHat = SIMD3(
        (quadrupoleMoment.0 * rHat).sum(),
        (quadrupoleMoment.1 * rHat).sum(),
        (quadrupoleMoment.2 * rHat).sum())
      let tensorProduct = (rHat * QrHat).sum()
      potential += 1 / 2 * tensorProduct / (rNorm * rNorm * rNorm)
    }
    estimates.append(potential)
    
    
    
    // Evaluate the potential directly.
    potential = .zero
    for charge in pointCharges {
      let position = unsafeBitCast(charge, to: SIMD3<Float>.self)
      let delta = samplePosition - position
      let deltaNorm = (delta * delta).sum().squareRoot()
      potential += charge.w / deltaNorm
    }
    estimates.append(potential)
    
    print(estimates)
    print(estimates.map { ($0 - estimates.last!).magnitude })
  }
}
