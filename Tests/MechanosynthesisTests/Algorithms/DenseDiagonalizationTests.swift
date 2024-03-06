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
}
