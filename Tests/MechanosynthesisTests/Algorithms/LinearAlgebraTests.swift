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
    
    // Predict what the QTAQ algorithm should produce.
    var expectedQ1 = [Float](repeating: 0, count: n * n)
    for columnID in 0..<n {
      // Find the length of the original vector.
      var norm: Float = .zero
      for rowID in 0..<n {
        let address = columnID * n + rowID
        let value = originalMatrixT[address]
        norm += value * value
      }
      norm = 1 / norm.squareRoot()
      
      // Overwrite the current value of Q1[columnID].
      for rowID in 0..<n {
        let address = rowID * n + columnID
        let value = originalMatrixT[address]
        expectedQ1[address] = value * norm
      }
      
      // Proceed with modified Gram-Schmidt.
      
    }
    
    print()
    print("expected Q1 (column-wise GS)")
    for rowID in 0..<n {
      for columnID in 0..<n {
        let address = rowID * n + columnID
        let value = expectedQ1[address]
        print(value, terminator: ", ")
      }
      print()
    }
    
    print()
    print("expected Q1^T (column-wise GS)")
    let expectedQ1T = Self.transpose(matrix: expectedQ1, n: n)
    for rowID in 0..<n {
      for columnID in 0..<n {
        let address = rowID * n + columnID
        let value = expectedQ1T[address]
        print(value, terminator: ", ")
      }
      print()
    }
    
    // Compare to a previous experiment, which did row-wise Gram-Schmidt.
    
    print()
    print("expected Q1 (row-wise GS)")
    print()
    
    /*
     original T
     4.0, -2.9999995, -4.768367e-08, 3.3378603e-07,
     -2.9999995, 3.3333328, -1.666667, -2.3841858e-07,
     -4.768367e-08, -1.666667, -1.32, 0.90666723,
     3.3378603e-07, -1.1920929e-07, 0.90666723, 1.9866666,

     expected Q1 (column-wise GS)
     0.8000001, -0.62705976, -2.0630448e-08, 1.5284792e-07,
     -0.59999996, 0.69673306, -0.7210872, -1.0917707e-07,
     -9.536735e-09, -0.34836665, -0.571101, 0.4151827,
     6.675721e-08, -2.4917119e-08, 0.39227164, 0.909738,

     expected Q1^T (column-wise GS)
     0.8000001, -0.59999996, -9.536735e-09, 6.675721e-08,
     -0.62705976, 0.69673306, -0.34836665, -2.4917119e-08,
     -2.0630448e-08, -0.7210872, -0.571101, 0.39227164,
     1.5284792e-07, -1.0917707e-07, 0.4151827, 0.909738,
     */
  }
}
