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
    var matrixA: [Float] = [
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
    
    // Display the matrix 'Q'.
    var Q1 = [Float](repeating: 0, count: n * n)
    for rowID in 0..<n {
      for columnID in 0..<n {
        let start: Float = (rowID == columnID) ? 1 : 0
        let rowPart = v1[rowID]
        let columnPart = v1[columnID]
        let value = start - 2 * rowPart * columnPart
        
        let address = rowID * n + columnID
        Q1[address] = value
      }
    }
    print()
    print("Q1 matrix:")
    let expectedQ1: [Float] = [
      1, 0, 0, 0,
      0, -1 / 3, 2 / 3, -2 / 3,
      0, 2 / 3, 2 / 3, 1 / 3,
      0, -2 / 3, 1 / 3, 2 / 3,
    ]
    for rowID in 0..<n {
      for columnID in 0..<n {
        let address = rowID * n + columnID
        let value = Q1[address]
        print(value, terminator: ", ")
        
        let expected = expectedQ1[address]
        XCTAssertEqual(value, expected, accuracy: 1e-6)
      }
      print()
    }
    
    // Multiply the matrices Q1 A Q1.
    // - Later, produce the same result through the more efficient O(n)
    //   transformations. Validate the output against this GEMM kernel.
    func sgemm(lhs: [Float], rhs: [Float]) -> [Float] {
      var output = [Float](repeating: 0, count: n * n)
      for rowID in 0..<n {
        for columnID in 0..<n {
          var dotProduct: Float = .zero
          for k in 0..<n {
            let addressLHS = rowID * n + k
            let addressRHS = k * n + columnID
            let valueLHS = lhs[addressLHS]
            let valueRHS = rhs[addressRHS]
            dotProduct += valueLHS * valueRHS
          }
          let address = rowID * n + columnID
          output[address] = dotProduct
        }
      }
      return output
    }
    
    let AQ = sgemm(lhs: matrixA, rhs: Q1)
    print()
    print("AQ matrix:")
    let expectedAQ: [Float] = [
      4, -3, 0, 0,
      1, -4 / 3, 5 / 3, -2 / 3,
      -2, 10 / 3, 4 / 3, -1 / 3,
      2, -1, -1, -2,
    ]
    for rowID in 0..<n {
      for columnID in 0..<n {
        let address = rowID * n + columnID
        let value = AQ[address]
        print(value, terminator: ", ")
        
        let expected = expectedAQ[address]
        XCTAssertEqual(value, expected, accuracy: 1e-6)
      }
      print()
    }
    
    let QAQ = sgemm(lhs: Q1, rhs: AQ)
    print()
    print("QAQ matrix:")
    var expectedQAQ: [Float] = [
      4, -3, 0, 0,
      -3, 10 / 3, 1, 4 / 3,
      0, 1, 5 / 3, -4 / 3,
      0, 4 / 3, -4 / 3, -1,
    ]
    for rowID in 0..<n {
      for columnID in 0..<n {
        let address = rowID * n + columnID
        let value = QAQ[address]
        print(value, terminator: ", ")
        
        let expected = expectedQAQ[address]
        XCTAssertEqual(value, expected, accuracy: 1e-6)
      }
      print()
    }
    
    // Create a standalone function for tridiagonalizing a matrix.
    func tridiagonalize(
      matrix originalMatrix: [Float],
      n: Int
    ) -> [Float] {
      /*
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
       */
      
      // This requires that n > 1.
      var currentMatrixA = originalMatrix
      for k in 0..<n - 2 {
        // Determine 'α' and 'r'.
        let address_k1k = (k + 1) * n + k
        let ak1k = currentMatrixA[address_k1k]
        var alpha: Float = .zero
        
        for rowID in (k + 1)..<n {
          let columnID = k
          let address = rowID * n + columnID
          let value = currentMatrixA[address]
          alpha += value * value
        }
        alpha.formSquareRoot()
        alpha *= (ak1k >= 0) ? -1 : 1
        
        var r = alpha * alpha - ak1k * alpha
        r = (r / 2).squareRoot()
        print(ak1k) // 1.0
        print(alpha) // -3.0
        print(r) // 2.449
        
        // Construct 'v'.
        var v = [Float](repeating: 0, count: n)
        v[k + 1] = (ak1k - alpha) / (2 * r)
        for vectorLane in (k + 2)..<n {
          let matrixAddress = vectorLane * n + k
          let matrixValue = currentMatrixA[matrixAddress]
          v[vectorLane] = matrixValue / (2 * r)
        }
        print(v) // [0.0, 0.81649655, -0.40824828, 0.40824828]
        
        // Operation 1: gemv(A, v)
        var operationResult1 = [Float](repeating: 0, count: n)
        for rowID in 0..<n {
          var dotProduct: Float = .zero
          for columnID in 0..<n {
            let matrixAddress = rowID * n + columnID
            let vectorAddress = columnID
            let matrixValue = currentMatrixA[matrixAddress]
            let vectorValue = v[vectorAddress]
            dotProduct += matrixValue * vectorValue
          }
          operationResult1[rowID] = dotProduct
        }
        
        // Operation 2: scatter(..., v)
        for rowID in 0..<n {
          for columnID in 0..<n {
            let rowValue = operationResult1[rowID]
            let columnValue = v[columnID]
            let matrixAddress = rowID * n + columnID
            currentMatrixA[matrixAddress] -= 2 * rowValue * columnValue
          }
        }
        
        // Check that the current matrix (AQ) matches expectations.
        print()
        print("AQ matrix:")
        for rowID in 0..<n {
          for columnID in 0..<n {
            let address = rowID * n + columnID
            let value = currentMatrixA[address]
            print(value, terminator: ", ")
          }
          print()
        }
        
        // Operation 3: gemv(vT, A)
        var operationResult3 = [Float](repeating: 0, count: n)
        for columnID in 0..<n {
          var dotProduct: Float = .zero
          for rowID in 0..<n {
            let vectorAddress = rowID
            let matrixAddress = rowID * n + columnID
            let vectorValue = v[vectorAddress]
            let matrixValue = currentMatrixA[matrixAddress]
            dotProduct += vectorValue * matrixValue
          }
          operationResult3[columnID] = dotProduct
        }
        
        // Operation 4: scatter(v, ...)
        for rowID in 0..<n {
          for columnID in 0..<n {
            let rowValue = v[rowID]
            let columnValue = operationResult3[columnID]
            let matrixAddress = rowID * n + columnID
            currentMatrixA[matrixAddress] -= 2 * rowValue * columnValue
          }
        }
        
        // Check that the current matrix (QAQ) matches expectations.
        print()
        print("QAQ matrix:")
        for rowID in 0..<n {
          for columnID in 0..<n {
            let address = rowID * n + columnID
            let value = currentMatrixA[address]
            print(value, terminator: ", ")
          }
          print()
        }
      }
      return []
    }
    
    print()
    print("tridiagonalizing")
    let hessenbergForm = tridiagonalize(matrix: matrixA, n: 4)
    
    print()
    print("hessenberg form:")
    print(hessenbergForm)
  }
}
