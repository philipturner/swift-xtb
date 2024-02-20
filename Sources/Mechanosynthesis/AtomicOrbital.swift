//
//  AtomicOrbital.swift
//
//
//  Created by Philip Turner on 2/12/24.
//

import Numerics

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
  var radialPart: (Float) -> Float
  var angularPart: (Float, Float, Float, Float) -> Float
  
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
      var output: Float = 1
      for _ in 0..<l {
        output *= shellRadiusPart
      }
      
      output *= Float.exp(-shellRadiusPart / 2)
      output *= L(shellRadiusPart)
      return output
    }
    angularPart = cubicHarmonic(l: l, m: m)
  }
  
  // Enter the position in a coordinate space where the nucleus is the origin.
  func waveFunction(position: SIMD3<Float>) -> Float {
    let r = (position * position).sum().squareRoot()
    let R = radialPart(r)
    let Y = angularPart(position.x, position.y, position.z, r)
    return normalizationFactor * R * Y
  }
}

func createShellOccupations(Z: Int) -> [Int] {
  // Output is zero-indexed, starting with the nonexistent zero shell. It spans
  // from n=0 to n=7.
  var shellOccupations = [Int](repeating: 0, count: 1 + 7)
  
  var cursorZ: Int = 0
  func subShell(n: Int, occupancy: Int) {
    if Z > cursorZ {
      shellOccupations[n] += min(Z - cursorZ, occupancy)
    }
    cursorZ += occupancy
  }
  
  // First period.
  subShell(n: 1, occupancy: 2)
  
  // Second period.
  subShell(n: 2, occupancy: 2)
  subShell(n: 2, occupancy: 6)
  
  // Third period.
  subShell(n: 3, occupancy: 2)
  subShell(n: 3, occupancy: 6)
  
  // Fourth period.
  subShell(n: 4, occupancy: 2)
  subShell(n: 3, occupancy: 10)
  subShell(n: 4, occupancy: 6)
  
  // Fifth period.
  subShell(n: 5, occupancy: 2)
  subShell(n: 4, occupancy: 10)
  subShell(n: 5, occupancy: 6)
  
  // Sixth period.
  subShell(n: 6, occupancy: 2)
  subShell(n: 4, occupancy: 14)
  subShell(n: 5, occupancy: 10)
  subShell(n: 6, occupancy: 6)
  
  // Seventh period.
  subShell(n: 7, occupancy: 2)
  subShell(n: 5, occupancy: 14)
  subShell(n: 6, occupancy: 10)
  subShell(n: 7, occupancy: 6)
  
  // Eighth period.
  if Z > 118 {
    fatalError("Eighth period elements are not supported.")
  }
  
  return shellOccupations
}

// Returns the effective charge for each electron shell.
//
// Gives ions the same effective charge as neutral, minimally spin-polarized
// atoms.
func createEffectiveCharges(Z: Int) -> [Float] {
  let occupations = createShellOccupations(Z: Z)
  var coreCharge = Z
  var effectiveCharges: [Float] = []
  for occupation in occupations {
    effectiveCharges.append(Float(coreCharge))
    coreCharge -= occupation
  }
  
  // Algorithm:
  // - break up the effective charge into different components
  // - fully shielded (0)
  // - partially shielded (sqrt(Z))
  // - not shielded (Z)
  for i in effectiveCharges.indices {
    var valenceCharge = Float(occupations[i])
    let coreCharge = effectiveCharges[i] - valenceCharge
    valenceCharge.formSquareRoot()
    effectiveCharges[i] = valenceCharge + coreCharge
  }
  
  return effectiveCharges
}
