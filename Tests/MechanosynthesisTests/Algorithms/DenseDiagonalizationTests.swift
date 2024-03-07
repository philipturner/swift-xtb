import XCTest
import Accelerate // Gate out once there are Swift kernels for CPUs without AMX.
import Numerics
import QuartzCore

// TODO: Create a custom direct diagonalization kernel. Start with the
// simplest one known, then move on to more advanced solvers.
final class DenseDiagonalizationTests: XCTestCase {
  // Reproduce the numerical procedure for Givens rotations.
  func testGivensRotation() {
    // MATLAB code from Wikipedia translated to optimal assembly code.
    // https://en.wikipedia.org/wiki/Givens_rotation#Stable_calculation
    func fastGivensRotation(a: Float, b: Float) -> (c: Float, s: Float) {
      let signA: Float = (a > 0) ? 1 : -1
      let signB: Float = (b > 0) ? 1 : -1
      let magnitudeA = a.magnitude
      let magnitudeB = b.magnitude
      
      var c: Float = 1
      var s: Float = 1
      if b == 0 {
        s = 0
      } else if a == 0 {
        c = 0
      } else {
        var t = min(magnitudeA, magnitudeB) / max(magnitudeA, magnitudeB)
        t *= (a * b > 0) ? 1 : -1
        c *= (magnitudeA > magnitudeB) ? 1 : t
        s *= (magnitudeA > magnitudeB) ? t : 1
        
        let u = 1 / (1 + t * t).squareRoot()
        c *= u
        s *= u
      }
      
      let multiplier = (magnitudeA > magnitudeB) ? signA : signB
      c *= multiplier
      s *= -multiplier
      return (c, s)
    }
    
    func checkGivensRotation(a: Float, b: Float) {
      let r = (a * a + b * b).squareRoot()
      let c = a / r
      let s = -b / r
      let (fastC, fastS) = fastGivensRotation(a: a, b: b)
      #if false
      print()
      print(c, "->", fastC)
      print(s, "->", fastS)
      print(r)
      #endif
      
      if r.isInfinite {
        let ratioBefore = a / b
        let ratioAfter = fastC / fastS
        
        if ratioBefore.isInfinite {
          XCTAssert(ratioAfter.isInfinite)
        } else if ratioBefore.isSubnormal || ratioBefore == 0 {
          XCTAssert(ratioAfter.isSubnormal || ratioAfter == 0)
        } else {
          XCTAssertEqual(ratioBefore / ratioAfter, -1, accuracy: 1e-6)
        }
      } else if r.isSubnormal || r == 0 {
        XCTAssert(c.isInfinite || c.isNaN)
        XCTAssert(s.isInfinite || s.isNaN)
      } else if r.isFinite {
        XCTAssertEqual(c, fastC, accuracy: 1e-6, "r=\(r)")
        XCTAssertEqual(s, fastS, accuracy: 1e-6, "r=\(r)")
      } else {
        fatalError("This should never happen.")
      }
    }
    
    func checkAbsoluteValuePair(first: Float, second: Float) {
      checkGivensRotation(a: first, b: second)
      checkGivensRotation(a: first, b: -second)
      checkGivensRotation(a: -first, b: second)
      checkGivensRotation(a: -first, b: -second)
    }
    
    // Test typical numbers.
    checkAbsoluteValuePair(first: 2.5237882239, second: 3.51613667762347)
    checkAbsoluteValuePair(first: 3.51613667762347, second: 2.5237882239)
    checkAbsoluteValuePair(first: 2.5237882239, second: 3.51613667762347e2)
    checkAbsoluteValuePair(first: 3.51613667762347e2, second: 2.5237882239)
    checkAbsoluteValuePair(first: 2.5237882239, second: 3.51613667762347e6)
    checkAbsoluteValuePair(first: 3.51613667762347e6, second: 2.5237882239)
    checkAbsoluteValuePair(first: 2.5237882239, second: 3.51613667762347e12)
    checkAbsoluteValuePair(first: 3.51613667762347e12, second: 2.5237882239)
    checkAbsoluteValuePair(first: 2.5, second: 3.5)
    checkAbsoluteValuePair(first: 3.5, second: 2.5)
    checkAbsoluteValuePair(first: 0.3, second: 0.7)
    checkAbsoluteValuePair(first: 0.7, second: 0.3)
    checkAbsoluteValuePair(first: 0.3, second: 0.0)
    checkAbsoluteValuePair(first: 0.0, second: 0.3)
    checkAbsoluteValuePair(first: 2.5, second: 0.0)
    checkAbsoluteValuePair(first: 0.0, second: 2.5)
    
    // Test numbers that overflow/underflow.
    checkAbsoluteValuePair(first: 2e25, second: 0.0)
    checkAbsoluteValuePair(first: 0.0, second: 2e25)
    checkAbsoluteValuePair(first: 2e25, second: 0.7)
    checkAbsoluteValuePair(first: 0.7, second: 2e25)
    checkAbsoluteValuePair(first: 0.3, second: 2e-25)
    checkAbsoluteValuePair(first: 2e-25, second: 0.3)
    checkAbsoluteValuePair(first: 0.0, second: 2e-25)
    checkAbsoluteValuePair(first: 2e-25, second: 0.0)
    checkAbsoluteValuePair(first: 1e36, second: 0)
    checkAbsoluteValuePair(first: 0, second: 1e36)
    checkAbsoluteValuePair(first: 1e36, second: 1e-36)
    checkAbsoluteValuePair(first: 1e-36, second: 1e36)
  }
  
  // Reproduce the example on Wikipedia for tridiagonalization through
  // Householder transformations.
  // https://en.wikipedia.org/wiki/Householder_transformation#Tridiagonalization
  func testTridiagonalization() throws {
    // The matrix A is implicitly stored in row-major form. However, it's
    // symmetric, so the ordering doesn't matter.
    let n: Int = 4
    var matrixA = [Float](repeating: 0, count: n * n)
    matrixA = [
      4, 1, -2, 2,
      1, 2, 0, 1,
      -2, 0, 3, -2,
      2, 1, -2, -1
    ]
    
    // Check that the original matrix is symmetric.
    for rowID in 0..<n {
      for columnID in rowID..<n {
        let upperAddress = rowID * n + columnID
        let lowerAddress = columnID * n + rowID
        let upperValue = matrixA[upperAddress]
        let lowerValue = matrixA[lowerAddress]
        XCTAssertEqual(upperValue, lowerValue)
      }
    }
    
    #if true
    // Matrix entries are 1-indexed on Wikipedia.
    let a21 = matrixA[1 * n + 0]
    var alpha: Float = .zero
    for j in 1..<n {
      let value = matrixA[j * n + 0]
      alpha += value * value
    }
    alpha.formSquareRoot()
    alpha *= (a21 >= 0) ? -1 : 1
    
    var r = alpha * alpha - a21 * alpha
    r = (r / 2).squareRoot()
    print(a21)
    print(alpha)
    print(r)
    
    // Construct the vector 'v'.
    var v1 = [Float](repeating: 0, count: n)
    v1[0] = 0
    v1[1] = (a21 - alpha) / (2 * r)
    for k in 2..<n {
      let value = matrixA[k * n + 0]
      v1[k] = value / (2 * r)
    }
    print(v1)
    
    /*
     1.0
     -3.0
     2.4494898
     [0.0, 0.81649655, -0.40824828, 0.40824828]

     Q1 matrix:
     1.0, 0.0, 0.0, 0.0,
     0.0, -0.33333325, 0.6666666, -0.6666666,
     0.0, 0.6666666, 0.6666667, 0.3333333,
     0.0, -0.6666666, 0.3333333, 0.6666667,
     */
    #else
    var xNorm: Float = .zero
    for rowID in 0..<n {
      let address = rowID * n + 0
      let xValue = matrixA[address]
      xNorm += xValue * xValue
    }
    xNorm.formSquareRoot()
    xNorm *= (matrixA[0] > 0) ? 1 : -1
    print("x-norm:", xNorm)
    
    var u = [Float](repeating: 0, count: n)
    for rowID in 0..<n {
      let address = rowID * n + 0
      let xValue = matrixA[address]
      let eValue: Float = (rowID == 0) ? 1 : 0
      u[rowID] = xValue - xNorm * eValue
    }
    print("u:", u)
    var uNorm: Float = .zero
    for rowID in 0..<n {
      let uValue = u[rowID]
      uNorm += uValue * uValue
    }
    uNorm.formSquareRoot()
    print("u-norm:", uNorm)
    
    var v1 = [Float](repeating: 0, count: n)
    for rowID in 0..<n {
      let uValue = u[rowID]
      v1[rowID] = uValue / uNorm
    }
    
    /*
     x-norm: 5.0
     u: [-1.0, 1.0, -2.0, 2.0]
     u-norm: 3.1622777

     Q1 matrix:
     0.8, 0.2, -0.4, 0.4,
     0.2, 0.8, 0.4, -0.4,
     -0.4, 0.4, 0.19999999, 0.8,
     0.4, -0.4, 0.8, 0.19999999,
     */
    #endif
    
    // Display the matrix 'Q'.
    var Q1 = [Float](repeating: 0, count: n * n)
    for rowID in 0..<n {
      for columnID in 0..<n {
        let start: Float = (rowID == columnID) ? 1 : 0
        let rowPart = v1[rowID]
        let columnPart = v1[columnID]
        let value = start - 2 * rowPart * columnPart
//        let value = rowPart * columnPart
        
        let address = rowID * n + columnID
        Q1[address] = value
      }
    }
    print()
    print("Q1 matrix:")
    for rowID in 0..<n {
      for columnID in 0..<n {
        let address = rowID * n + columnID
        let value = Q1[address]
        print(value, terminator: ", ")
      }
      print()
    }
  }
}
