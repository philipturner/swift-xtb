//
//  BasisFunction+Math.swift
//
//
//  Created by Philip Turner on 2/13/24.
//

func laguerrePolynomial(
  alpha: Float, n: Int
) -> (_ x: SIMD8<Float>) -> SIMD8<Float> {
  if n == 0 {
    return { _ in .one }
  } else if n > 0 {
    return { x in
      var secondLast: SIMD8<Float> = .one
      var last: SIMD8<Float> = 1 + alpha - x
      
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
) -> (
  _ x: SIMD8<Float>,
  _ y: SIMD8<Float>,
  _ z: SIMD8<Float>,
  _ r: SIMD8<Float>
) -> SIMD8<Float> {
  var factorial: Int = 1
  for i in 0...l {
    factorial *= (i * 2 + 1)
  }
  let Nc = (Float(factorial) / (4 * Float.pi)).squareRoot()
  
  // The basis set from Wikipedia, which contains up to f-orbitals:
  // https://en.wikipedia.org/wiki/Cubic_harmonic
  if l == 0 {
    return { _, _, _, _ in
      var output = Nc / SIMD8<Float>.one
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
      case 0:
        output *= 3 * z * z - r * r
        output /= 2 * Float(3).squareRoot()
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
      
      // Prevent the compiler from taking a long time to type-check this.
      let x2: SIMD8<Float> = x * x
      let y2: SIMD8<Float> = y * y
      let z2: SIMD8<Float> = z * z
      
      switch m {
      case 0:
        let expr = 2 * z2 - 3 * x2 - 3 * y2
        output *= z * expr
        output /= 2 * Float(15).squareRoot()
      case -1:
        output *= x * (4 * z2 - x2 - y2)
        output /= 2 * Float(10).squareRoot()
      case 1:
        output *= y * (4 * z2 - x2 - y2)
        output /= 2 * Float(10).squareRoot()
      case -2:
        output *= x * y * z
      case 2:
        output *= z * (x2 - y2) / 2
      case -3:
        output *= x * (x2 - 3 * y2)
        output /= 2 * Float(6).squareRoot()
      case 3:
        output *= y * (3 * x2 - y2)
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
