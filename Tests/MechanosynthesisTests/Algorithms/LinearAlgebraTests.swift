import XCTest
import Accelerate // Gate out once there are Swift kernels for CPUs without AMX.
import Numerics
import QuartzCore

final class LinearAlgebraTests: XCTestCase {
  // MARK: - Linear Algebra Functions
  
  // Multiplies two square matrices.
  // - Mixed precision instead of single precision.
  static func matrixMultiply(
    matrixA: [Float], transposeA: Bool = false,
    matrixB: [Float], transposeB: Bool = false,
    n: Int
  ) -> [Float] {
    var matrixC = [Float](repeating: 0, count: n * n)
    for rowID in 0..<n {
      for columnID in 0..<n {
        var dotProduct: Double = .zero
        for k in 0..<n {
          var value1: Float
          var value2: Float
          if !transposeA {
            value1 = matrixA[rowID * n + k]
          } else {
            value1 = matrixA[k * n + rowID]
          }
          if !transposeB {
            value2 = matrixB[k * n + columnID]
          } else {
            value2 = matrixB[columnID * n + k]
          }
          dotProduct += Double(value1 * value2)
        }
        matrixC[rowID * n + columnID] = Float(dotProduct)
      }
    }
    return matrixC
  }
  
  // Forms an orthogonal basis of the matrix's rows.
  // - Mixed precision instead of single precision.
  // - Classical GS instead of modified GS.
  static func rowGramSchmidt(
    matrix originalMatrix: [Float], n: Int
  ) -> [Float] {
    // Operate on the (output) matrix in-place.
    var matrix = originalMatrix
    
    func normalize(electronID: Int) {
      var norm: Double = .zero
      for cellID in 0..<n {
        let value = matrix[electronID * n + cellID]
        norm += Double(value * value)
      }
      
      let normalizationFactor = 1 / norm.squareRoot()
      for cellID in 0..<n {
        var value = Double(matrix[electronID * n + cellID])
        value *= normalizationFactor
        matrix[electronID * n + cellID] = Float(value)
      }
    }
    
    for electronID in 0..<n {
      // Normalize the vectors before taking dot products.
      normalize(electronID: electronID)
    }
    
    for electronID in 0..<n {
      // Determine the magnitude of components parallel to previous vectors.
      var dotProducts = [Double](repeating: .zero, count: n)
      for neighborID in 0..<electronID {
        var dotProduct: Double = .zero
        for cellID in 0..<n {
          let value1 = matrix[electronID * n + cellID]
          let value2 = matrix[neighborID * n + cellID]
          dotProduct += Double(value1) * Double(value2)
        }
        dotProducts[neighborID] = dotProduct
      }
      
      // Subtract all components parallel to previous vectors.
      for cellID in 0..<n {
        var value1 = Double(matrix[electronID * n + cellID])
        for neighborID in 0..<electronID {
          let dotProduct = dotProducts[neighborID]
          let value2 = matrix[neighborID * n + cellID]
          value1 -= dotProduct * Double(value2)
        }
        matrix[electronID * n + cellID] = Float(value1)
      }
      
      // Rescale the orthogonal component to unit vector length.
      normalize(electronID: electronID)
    }
    
    return matrix
  }
  
  // Standalone function for tridiagonalizing a matrix.
  // - Performs several loop iterations with O(n^2) complexity each.
  // - Overall, the procedure has O(n^3) complexity.
  static func tridiagonalize(
    matrix originalMatrix: [Float],
    n: Int
  ) -> [Float] {
    var currentMatrixA = originalMatrix
    
    // This requires that n > 1.
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
      
      // Construct 'v'.
      var v = [Float](repeating: 0, count: n)
      v[k + 1] = (ak1k - alpha) / (2 * r)
      for vectorLane in (k + 2)..<n {
        let matrixAddress = vectorLane * n + k
        let matrixValue = currentMatrixA[matrixAddress]
        v[vectorLane] = matrixValue / (2 * r)
      }
      
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
    }
    return currentMatrixA
  }
  
  // Returns the transpose of a square matrix.
  static func transpose(matrix: [Float], n: Int) -> [Float] {
    var output = [Float](repeating: 0, count: n * n)
    for rowID in 0..<n {
      for columnID in 0..<n {
        let oldAddress = columnID * n + rowID
        let newAddress = rowID * n + columnID
        output[newAddress] = matrix[oldAddress]
      }
    }
    return output
  }
  
  // MARK: - Tests
  
  // Test the QR algorithm for finding eigenvalues of a tridiagonal matrix.
  func testQRAlgorithm() throws {
    let n: Int = 4
    let originalMatrixA: [Float] = [
      4, 1, -2, 2,
      1, 2, 0, 1,
      -2, 0, 3, -2,
      2, 1, -2, -1
    ]
    let originalMatrixT = Self.tridiagonalize(matrix: originalMatrixA, n: n)
    print()
    print("original T")
    for rowID in 0..<n {
      for columnID in 0..<n {
        let address = rowID * n + columnID
        let value = originalMatrixT[address]
        print(value, terminator: ", ")
      }
      print()
    }
    
    // Predict what the Q^T A Q / Q^T T Q algorithm should produce.
    var expectedQ1 = originalMatrixT
    var expectedQ1T = Self.transpose(matrix: expectedQ1, n: n)
    expectedQ1T = Self.rowGramSchmidt(matrix: expectedQ1T, n: n)
    expectedQ1 = Self.transpose(matrix: expectedQ1T, n: n)
    print()
    print("expected Q1")
    for rowID in 0..<n {
      for columnID in 0..<n {
        let address = rowID * n + columnID
        let value = expectedQ1[address]
        print(value, terminator: ", ")
      }
      print()
    }
    
    var expectedT2 = Self.matrixMultiply(
      matrixA: originalMatrixT, matrixB: expectedQ1, n: n)
    expectedT2 = Self.matrixMultiply(
      matrixA: expectedQ1T, matrixB: expectedT2, n: n)
    print()
    print("expected T2")
    for rowID in 0..<n {
      for columnID in 0..<n {
        let address = rowID * n + columnID
        let value = expectedT2[address]
        print(value, terminator: ", ")
      }
      print()
    }
    
    // Store the tridiagonal matrix in a compact form.
    let currentMatrixT = originalMatrixT
    var a = [Float](repeating: 0, count: n + 1)
    var bSquared = [Float](repeating: 0, count: n)
    for diagonalID in 0..<n {
      let matrixAddress = diagonalID * n + diagonalID
      let vectorAddress = diagonalID
      a[vectorAddress] = currentMatrixT[matrixAddress]
    }
    for subDiagonalID in 0..<n - 1 {
      let rowID = subDiagonalID
      let columnID = subDiagonalID + 1
      let matrixAddress = rowID * n + columnID
      let vectorAddress = rowID
      
      let bValue = currentMatrixT[matrixAddress]
      bSquared[vectorAddress] = bValue * bValue
    }
    
    // Perform one iteration of QR decomposition.
    // https://doi.org/10.1093/comjnl/6.1.99
    var converged = false
    for _ in 0..<10 {
      // Shift the matrix: T -> T - λI
      let λ = Self.shiftJZ(a: a, bSquared: bSquared, n: n)
      for diagonalID in 0..<n {
        a[diagonalID] -= λ
      }
      
      // Factorize Q, R and compute RQ.
      var u: Float = .zero
      var sSquared: Float = .zero
      var oneMinusSSquared: Float = 1
      for i in 0..<n {
        let γ = a[i] - u
        let pSquared = (γ * γ) / oneMinusSSquared
        let p2b2 = pSquared + bSquared[i]
        if i > 0 {
          bSquared[i - 1] = sSquared * (p2b2)
        }
        
        // Skip the loop iteration when a[i] becomes zero. This means we
        // have converged.
        if p2b2.magnitude < .leastNormalMagnitude {
          guard i == n - 1 else {
            fatalError("a[i] was zero, but i was not n - 1.")
          }
          print("About to be zero. Skipping iteration \(i).")
          converged = true
          continue
        }
        
        // Store 1 - s^2 in a way that's numerically stable.
        sSquared = bSquared[i] / (p2b2)
        oneMinusSSquared = pSquared / (p2b2)
        u = sSquared * (γ + a[i + 1])
        a[i] = γ + u
      }
      
      // Un-shift the matrix: RQ -> RQ + λI
      for diagonalID in 0..<n {
        a[diagonalID] += λ
      }
      
      print()
      print("actual a")
      for slotID in a.indices {
        print(slotID, "| a =", a[slotID])
      }
      
      print()
      print("actual b")
      for slotID in bSquared.indices {
        let value = bSquared[slotID]
        print(slotID, "| b^2 =", value, terminator: ", ")
        print("|b| =", value.squareRoot())
      }
    }
    guard converged else {
      fatalError("QR decomposition failed to converge after 10 trials.")
    }
  }
  
  // Compute the JZ shift for a QR iteration.
  //
  // Arguments:
  // - a (arbitrary length array)
  // - bSquared (arbitrary length array)
  // - n (integer)
  // Returns:
  // - λ (floating-point number)
  static func shiftJZ(a: [Float], bSquared: [Float], n: Int) -> Float {
    guard n - 2 >= 0 else {
      fatalError("Could not perform Wilkinson shift.")
    }
    let an = a[n - 1]
    let an1 = a[n - 2]
    let bSquaredn1 = bSquared[n - 2]
    
    // Solve the quadratic equation:
    // (a_{n - 1} - λ)(a_{n} - λ) = (b_{n - 1})^2
    let a: Float = 1
    let b: Float = -(an + an1)
    let c: Float = an * an1 - bSquaredn1
    let determinant = b * b - 4 * a * c
    guard determinant > 0 else {
      print(
        "Could not solve quadratic equation. Determinant was \(determinant).")
      return an
    }
    var solution1 = (-b + determinant.squareRoot()) / (2 * a)
    var solution2 = (-b - determinant.squareRoot()) / (2 * a)
    
    // Select the closest eigenvalue:
    // |a_{n} - λ| ≤ |b_{n - 1}| ≤ |a_{n - 1} - λ|
    func inequalityResidual(λ: Float) -> Float {
      let valueLeft = (an - λ).magnitude
      let valueMiddle = bSquaredn1.squareRoot()
      let valueRight = (an1 - λ).magnitude
      var inequalityResidual = max(0, valueLeft - valueMiddle)
      inequalityResidual += max(0, valueMiddle - valueRight)
      return inequalityResidual
    }
    let residual1 = inequalityResidual(λ: solution1)
    let residual2 = inequalityResidual(λ: solution2)
    
    // Store the Wilkinson shift in 'solution1'.
    if residual2 < residual1 {
      swap(&solution1, &solution2)
    }
    
    // Check whether the closest eigenvalue actually satisfies the
    // inequality above.
    let closenessResidual1 = (solution1 - an).magnitude
    let closenessResidual2 = (solution2 - an).magnitude
    guard closenessResidual1 < closenessResidual2 else {
      fatalError("Incorrect solution ordering.")
    }
    
    // Print the values obtained from the Rayleigh and Wilkinson shifts.
    // Return the Rayleigh shift for now.
    print("solution 1:", solution1)
    print("solution 2:", solution2)
    print("Rayleigh shift:", an)
    
    guard n - 3 >= 0 else {
      fatalError("Could not fetch b_{n - 2}.")
    }
    let bSquaredn2 = bSquared[n - 3]
    if bSquaredn2.squareRoot() >= (2 * bSquaredn1).squareRoot() {
      print("Returning Rayleigh shift.")
      return an
    } else {
      print("Returning Wilkinson shift.")
      return solution1
    }
  }
}
