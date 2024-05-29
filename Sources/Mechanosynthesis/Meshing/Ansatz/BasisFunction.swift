//
//  BasisFunction.swift
//
//
//  Created by Philip Turner on 2/12/24.
//

import Numerics

/// A configuration for a basis function.
public struct BasisFunctionDescriptor {
  /// Required. Effective nuclear charge.
  public var Z: Float?
  
  /// Required. Primary quantum number.
  public var n: Int?
  
  /// Required. Angular momentum quantum number.
  public var l: Int?
  
  /// Required. Magnetic quantum number.
  public var m: Int?
}

/// A basis function for a hydrogenic orbital.
public struct BasisFunction {
  /// Effective nuclear charge.
  public var Z: Float
  
  /// Primary quantum number.
  public var n: Int
  
  /// Angular momentum quantum number.
  public var l: Int
  
  /// Magnetic quantum number.
  public var m: Int
  
  @usableFromInline
  var normalizationFactor: Float
  
  @usableFromInline
  var L: (SIMD8<Float>) -> SIMD8<Float>
  
  @usableFromInline
  var Y: (
    SIMD8<Float>, SIMD8<Float>, SIMD8<Float>) -> SIMD8<Float>
  
  init(descriptor: BasisFunctionDescriptor) {
    guard let Z = descriptor.Z,
          let n = descriptor.n,
          let l = descriptor.l,
          let m = descriptor.m else {
      fatalError("Descriptor was incomplete.")
    }
    self.Z = Z
    self.n = n
    self.l = l
    self.m = m
    
    // Correction to achieve the same scaling behavior as actual atomic
    // orbitals: remove the division by 'n' from the shell part.
    let shellPart = Float(2 * Z)
    normalizationFactor = Float(factorial(n - l - 1))
    normalizationFactor /= Float(2 * n * factorial(n + l))
    normalizationFactor *= shellPart * shellPart * shellPart
    normalizationFactor.formSquareRoot()
    
    L = laguerrePolynomial(
      alpha: Float(2 * l + 1), n: n - l - 1)
    Y = cubicHarmonic(l: l, m: m)
  }
  
  /// Vectorized function for calculating the radial part.
  ///
  /// - Parameter r: The distance from the nucleus.
  @_transparent
  public func radialPart(r: SIMD8<Float>) -> SIMD8<Float> {
    var output: SIMD8<Float> = .init(repeating: normalizationFactor)
    
    // Define the shell-radius part.
    let shellPart = Float(2 * Z)
    let shellRadiusPart = shellPart * r
    for _ in 0..<l {
      output *= shellRadiusPart
    }
    
    // Hard-coded the logarithm of e in base 2.
    let input = -shellRadiusPart / 2 * 1.4426950408889607
    var expValue: SIMD8<Float> = .zero
    for laneID in 0..<8 {
      // Unsure how to vectorize transcendentals on all platforms.
      expValue[laneID] = Float.exp2(input[laneID])
    }
    output *= expValue
    output *= L(shellRadiusPart)
    return output
  }
  
  /// Vectorized function for calculating the angular part.
  ///
  /// - Parameter x: The x-coordinate of the vector separating the point from
  ///                the nucleus's position.
  /// - Parameter y: The y-coordinate of the vector separating the point from
  ///                the nucleus's position.
  /// - Parameter z: The z-coordinate of the vector separating the point from
  ///                the nucleus's position.
  @_transparent
  public func angularPart(
    x: SIMD8<Float>,
    y: SIMD8<Float>,
    z: SIMD8<Float>
  ) -> SIMD8<Float> {
    return Y(x, y, z)
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
