//
//  AtomicOrbital.swift
//
//
//  Created by Philip Turner on 2/12/24.
//

import Numerics
#if canImport(simd)
import simd
#endif

struct AtomicOrbitalDescriptor {
  // Effective nuclear charge.
  var Z: Float?
  
  // Primary quantum number.
  var n: Int?
  
  // Angular momentum quantum number.
  var l: Int?
  
  // Magnetic quantum number.
  var m: Int?
}

struct AtomicOrbital {
  var normalizationFactor: Float
  var radialPart: (SIMD8<Float>) -> SIMD8<Float>
  var angularPart: (
    SIMD8<Float>, SIMD8<Float>, SIMD8<Float>, SIMD8<Float>) -> SIMD8<Float>
  
  init(descriptor: AtomicOrbitalDescriptor) {
    guard let Z = descriptor.Z,
          let n = descriptor.n,
          let l = descriptor.l,
          let m = descriptor.m else {
      fatalError("Atomic orbital descriptor was invalid.")
    }
    
    // Correction to achieve the same scaling behavior as actual atomic
    // orbitals: remove the division by 'n' from the shell part.
    let shellPart = Float(2 * Z)
    normalizationFactor = Float(factorial(n - l - 1))
    normalizationFactor /= Float(2 * n * factorial(n + l))
    normalizationFactor *= shellPart * shellPart * shellPart
    normalizationFactor.formSquareRoot()
    
    let L = laguerrePolynomial(alpha: Float(2 * l + 1), n: n - l - 1)
    radialPart = { r in
      let shellRadiusPart = shellPart * r
      var output: SIMD8<Float> = .one
      for _ in 0..<l {
        output *= shellRadiusPart
      }
      
      // Unsure how to make transcendentals faster on other platforms yet.
      // It may become irrelevant when everything is ported to GPU.
      let input = -shellRadiusPart / 2 * 1.4426950408889607
      var expValue: SIMD8<Float> = .zero
      #if canImport(simd)
      expValue = exp2(input)
      #else
      for laneID in 0..<8 {
        expValue[laneID] = Float.exp2(input[laneID])
      }
      #endif
      output *= expValue
      output *= L(shellRadiusPart)
      return output
    }
    angularPart = cubicHarmonic(l: l, m: m)
  }
  
  // Enter the position in a coordinate space where the nucleus is the origin.
  func waveFunction(
    x: SIMD8<Float>,
    y: SIMD8<Float>,
    z: SIMD8<Float>
  ) -> SIMD8<Float> {
    let r = (x * x + y * y + z * z).squareRoot()
    let R = radialPart(r)
    let Y = angularPart(x, y, z, r)
    return normalizationFactor * R * Y
  }
}

// Processes the electron count in a single spin channel.
func createShellOccupations(electronCount: Int) -> [Int] {
  // The output is zero-indexed, starting with the nonexistent zero shell.
  var shellOccupations = [Int](repeating: 0, count: 8)
  
  var cursor: Int = 0
  func subShell(n: Int, occupancy: Int) {
    if electronCount > cursor {
      shellOccupations[n] += min(electronCount - cursor, occupancy)
    }
    cursor += occupancy
  }
  
  // First period.
  subShell(n: 1, occupancy: 1)
  
  // Second period.
  subShell(n: 2, occupancy: 1)
  subShell(n: 2, occupancy: 3)
  
  // Third period.
  subShell(n: 3, occupancy: 1)
  subShell(n: 3, occupancy: 3)
  
  // Fourth period.
  subShell(n: 4, occupancy: 1)
  subShell(n: 3, occupancy: 5)
  subShell(n: 4, occupancy: 3)
  
  // Fifth period.
  subShell(n: 5, occupancy: 1)
  subShell(n: 4, occupancy: 5)
  subShell(n: 5, occupancy: 3)
  
  // Sixth period.
  subShell(n: 6, occupancy: 1)
  subShell(n: 4, occupancy: 7)
  subShell(n: 5, occupancy: 5)
  subShell(n: 6, occupancy: 3)
  
  // Seventh period.
  subShell(n: 7, occupancy: 1)
  subShell(n: 5, occupancy: 7)
  subShell(n: 6, occupancy: 5)
  subShell(n: 7, occupancy: 3)
  
  // Eighth period.
  if cursor < electronCount {
    fatalError("Too many electrons to create a shell occupation.")
  }
  
  return shellOccupations
}

// Returns the effective charge for each electron shell.
func createEffectiveCharges(
  Z: UInt8, spinDownOccupations: [Int], spinUpOccupations: [Int]
) -> [Float] {
  // Break up the effective charge into different components:
  // - not shielded (Z)
  // - partially shielded (sqrt(Z))
  // - fully shielded (0)
  var notShieldedCharge = Int(Z)
  var effectiveCharges: [Float] = []
  for shellID in 0..<8 {
    let spinDownOccupation = spinDownOccupations[shellID]
    let spinUpOccupation = spinUpOccupations[shellID]
    let occupation = spinDownOccupation + spinUpOccupation
    notShieldedCharge -= occupation
    
    // Traveling down a period, the net charge in a particular shell
    // increases. This shrinks the atomic radius. The rate of shrinking
    // scales with the square root of column number.
    let partiallyShieldedCharge = Float(occupation).squareRoot()
    let charge = Float(notShieldedCharge) + partiallyShieldedCharge
    
    // Clamp the effective charge to at least +1.
    effectiveCharges.append(max(1, charge))
  }
  
  return effectiveCharges
}
