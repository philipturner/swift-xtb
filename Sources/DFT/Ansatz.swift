//
//  Ansatz.swift
//
//
//  Created by Philip Turner on 2/4/24.
//

import OpenCL

// For reference: an implementation of the ansatz that reproduces atomic
// radii very well.
// - We'll need to slightly modify it for core electrons. It only handles
//   the valence electrons for now, to reproduce experimental data.
// - Make the bottom term scale with Z_eff = Z_inner + sqrt(Z_valence).
// - Embed the CPU code as a reference implementation for the unit tests.
//   Assert that the GPU produces the same results.

#if false

func createGeometry() -> [Entity] {
  //  for Z in 1...118 {
  //    var n: Int
  //    var Zeff: Int
  //    if Z <= 2 {
  //      n = 1
  //      Zeff = Z
  //    } else if Z <= 10 {
  //      n = 2
  //      Zeff = Z - 2
  //    } else if Z <= 18 {
  //      n = 3
  //      Zeff = Z - 10
  //    } else if Z <= 36 {
  //      n = 4
  //      if Z <= 20 {
  //        Zeff = Z - 18
  //      } else if Z <= 30 {
  //        Zeff = 2
  //      } else {
  //        Zeff = 2 + (Z - 30)
  //      }
  //    } else if Z <= 54 {
  //      n = 5
  //      if Z <= 38 {
  //        Zeff = Z - 36
  //      } else if Z <= 48 {
  //        Zeff = 2
  //      } else {
  //        Zeff = 2 + (Z - 48)
  //      }
  //    } else if Z <= 86 {
  //      n = 6
  //      if Z <= 56 {
  //        Zeff = Z - 54
  //      } else if Z <= 80 {
  //        Zeff = 2
  //      } else {
  //        Zeff = 2 + (Z - 80)
  //      }
  //    } else {
  //      n = 7
  //      if Z <= 88 {
  //        Zeff = Z - 86
  //      } else if Z <= 112 {
  //        Zeff = 2
  //      } else {
  //        Zeff = 2 + (Z - 112)
  //      }
  //    }
  //
  //    var correction = Float(n) * Float(n) / Float(Zeff)
  //    correction = correction.squareRoot()
  //    var normalization = correction * correction * correction
  //    normalization = normalization.squareRoot()
  
  let n: Int = 4
  let Zeff: Int = 1
  
  let numbers = QuantumNumbers(n: n, l: 3, m: 3, Z: Zeff)
  var integral: Double = .zero
  let maxRadius: Double = 50
  let dr: Double = 0.25
  
  // TODO: Create a rule for positioning any electron within an atom core.
  // Graph the trend in 's' electron radii across the different 'Z' and 'n'
  // quantum numbers. Next, analyze the overlap matrix of the 's' electrons.
  // After ensuring that every element looks sound, compute the radii and
  // overlap matrix of the all-electron Au atom.
  
  var z: Double = -maxRadius
  while z <= maxRadius {
    defer { z += dr }
    
    var y: Double = -maxRadius
    while y <= maxRadius {
      defer { y += dr }
      
      var x: Double = -maxRadius
      while x < maxRadius {
        defer { x += dr }
        
        let position = SIMD3(Float(x), Float(y), Float(z))
        let waveFunction = Double(hydrogenWaveFunction(
          numbers: numbers, position: position))
        
        let r = Double((position * position).sum().squareRoot())
        _ = r
        let observable: Double = waveFunction * z * z * waveFunction
        if observable.isNaN {
          print("NaN at \(x), \(y), \(z)")
          continue
        }
        integral += observable * dr * dr * dr
      }
    }
  }
  
  print("\(integral.squareRoot()),")
  //  }
  exit(0)
}

// MARK: - Hydrogen Wave Function

func laguerrePolynomial(
  alpha: Float, n: Int
) -> (_ x: Float) -> Float {
  if n == 0 {
    return { _ in 1 }
  } else if n > 0 {
    return { x in
      var secondLast: Float = 1
      var last: Float = 1 + alpha - x
      
      for k in 1..<n {
        let coeffLeft = Float(2 * k + 1) + alpha - x
        let coeffRight = -(Float(k) + alpha)
        let numerator = coeffLeft * last + coeffRight * secondLast
        let denominator = Float(k + 1)
        secondLast = last
        last = numerator / denominator
      }
      return last
    }
  }
  
  fatalError("Unsupported value for n.")
}

func cubicHarmonic(
  l: Int, m: Int
) -> (_ x: Float, _ y: Float, _ z: Float, _ r: Float) -> Float {
  var factorial: Int = 1
  for i in 0...l {
    factorial *= (i * 2 + 1)
  }
  let Nc = (Float(factorial) / (4 * Float.pi)).squareRoot()
  
  if l == 0 {
    return { _, _, _, _ in
      var output = Nc
      switch m {
      case 0: output *= 1
      default: fatalError("Invalid value for m.")
      }
      return output
    }
  } else if l == 1 {
    return { x, y, z, r in
      var output = Nc / r
      switch m {
      case 0: output *= z
      case -1: output *= x
      case 1: output *= y
      default: fatalError("Invalid value for m.")
      }
      return output
    }
  } else if l == 2 {
    return { x, y, z, r in
      var output = Nc / (r * r)
      switch m {
      case 0: output *= (3 * z * z - r * r) / (2 * Float(3).squareRoot())
      case -1: output *= x * z
      case 1: output *= y * z
      case -2: output *= x * y
      case 2: output *= (x * x - y * y) / 2
      default: fatalError("Invalid value for m.")
      }
      return output
    }
  } else if l == 3 {
    return { x, y, z, r in
      var output = Nc / (r * r * r)
      switch m {
      case 0:
        output *= z * (2 * z * z - 3 * x * x - 3 * y * y)
        output /= 2 * Float(15).squareRoot()
      case -1:
        output *= x * (4 * z * z - x * x - y * y)
        output /= 2 * Float(10).squareRoot()
      case 1:
        output *= y * (4 * z * z - x * x - y * y)
        output /= 2 * Float(10).squareRoot()
      case -2:
        output *= x * y * z
      case 2:
        output *= z * (x * x - y * y) / 2
      case -3:
        output *= x * (x * x - 3 * y * y)
        output /= 2 * Float(6).squareRoot()
      case 3:
        output *= y * (3 * x * x - y * y)
        output /= 2 * Float(6).squareRoot()
      default:
        fatalError("Invalid value for m.")
      }
      return output
    }
  }
  
  fatalError("Unsupported value for l.")
}

func factorial(_ x: Int) -> Int {
  guard x >= 0 else {
    fatalError("Cannot take factorial of negative number.")
  }
  if x == 0 {
    return 1
  } else {
    var output = x
    var counter = x - 1
    while counter > 0 {
      output *= counter
      counter -= 1
    }
    return output
  }
}

struct QuantumNumbers {
  var n: Int
  var l: Int
  var m: Int
  var Z: Int
}

func hydrogenWaveFunction(
  numbers: QuantumNumbers,
  position: SIMD3<Float>
) -> Float {
  let R = { (r: Float) -> Float in
    let numerator = factorial(numbers.n - numbers.l - 1)
    let denominator = 2 * numbers.n * factorial(numbers.n + numbers.l)
    var normalizationFactor = Float(numerator) / Float(denominator)
    
    let shellPart = Float(2 * numbers.Z) / Float(numbers.n)
    normalizationFactor *= shellPart * shellPart * shellPart
    normalizationFactor.formSquareRoot()
    
    let shellRadiusPart = shellPart * r
    let L = laguerrePolynomial(
      alpha: Float(2 * numbers.l + 1),
      n: numbers.n - numbers.l - 1)
    
    return normalizationFactor
    * exp(-shellRadiusPart / 2)
    * pow(shellRadiusPart, Float(numbers.l))
    * L(shellRadiusPart)
  }
  
  let r = (position * position).sum().squareRoot()
  let Y = cubicHarmonic(l: numbers.l, m: numbers.m)
  let magnitude = R(r) * Y(position.x, position.y, position.z, r)
  let parity = pow(-1, Float(numbers.l))
  return parity * magnitude
}

#endif
