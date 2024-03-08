import XCTest
import Accelerate
import Numerics
import QuartzCore

final class LinearAlgebraTests: XCTestCase {
  // MARK: - Linear Algebra Functions
  
  // Multiplies two square matrices.
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
  
  // Forms an orthogonal basis of the matrix's columns.
  static func modifiedGramSchmidt(
    matrix originalMatrix: [Float], n: Int
  ) -> [Float] {
    // Operate on the output matrix in-place.
    var matrix = originalMatrix
    
    func normalize(electronID: Int) {
      var norm: Double = .zero
      for cellID in 0..<n {
        let value = matrix[cellID * n + electronID]
        norm += Double(value * value)
      }
      
      let normalizationFactor = 1 / norm.squareRoot()
      for cellID in 0..<n {
        var value = matrix[cellID * n + electronID]
        value *= Float(normalizationFactor)
        matrix[cellID * n + electronID] = value
      }
    }
    
    for electronID in 0..<n {
      // Normalize the vectors before taking dot products.
      normalize(electronID: electronID)
    }
    
    for electronID in 0..<n {
      for neighborID in 0..<electronID {
        // Determine the magnitude of the parallel component.
        var dotProduct: Double = .zero
        for cellID in 0..<n {
          let value1 = matrix[cellID * n + electronID]
          let value2 = matrix[cellID * n + neighborID]
          dotProduct += Double(value1) * Double(value2)
        }
        
        // Subtract the parallel component.
        for cellID in 0..<n {
          var value1 = matrix[cellID * n + electronID]
          let value2 = matrix[cellID * n + neighborID]
          value1 -= Float(dotProduct) * value2
          matrix[cellID * n + electronID] = value1
        }
      }
      
      // Rescale the orthogonal component to unit vector length.
      normalize(electronID: electronID)
    }
    
    return matrix
  }
  
  // Reduce the matrix to tridiagonal form.
  // - TODO: Return the Householder transformations for back-transforming the
  //   eigenvectors.
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
  
  // Use the utility function from LAPACK to diagonalize a tridiagonal matrix.
  func testDivideAndConquer() throws {
    // TODO: Also test various 7x7 matrices from previous calculations,
    // especially the ones with degenerate eigenspaces.
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
    
    // Store the tridiagonal matrix in a compact form.
    var D = [Float](repeating: 0, count: n)
    var E = [Float](repeating: 0, count: n - 1)
    for diagonalID in 0..<n {
      let matrixAddress = diagonalID * n + diagonalID
      let vectorAddress = diagonalID
      D[vectorAddress] = originalMatrixT[matrixAddress]
    }
    for subDiagonalID in 0..<n - 1 {
      let rowID = subDiagonalID
      let columnID = subDiagonalID + 1
      let matrixAddress = rowID * n + columnID
      let vectorAddress = rowID
      E[vectorAddress] = originalMatrixT[matrixAddress]
    }
    
    // Query the workspace size.
    var JOBZ = CChar(Character("V").asciiValue!)
    var N = Int32(n)
    var LDZ = Int32(n)
    var WORK = [Float](repeating: 0, count: 1)
    var LWORK = Int32(-1)
    var IWORK = [Int32](repeating: 0, count: 1)
    var LIWORK = Int32(-1)
    var INFO = Int32(0)
    sstevd_(
      &JOBZ, &N, nil, nil, nil, &LDZ, &WORK, &LWORK, &IWORK, &LIWORK, &INFO)
    print("INFO:", INFO)
    print("LWORK:", WORK[0])
    print("LIWORK:", IWORK[0])
    print("D:", D)
    print("E:", E)
    
    // Call into LAPACK.
    var Z = [Float](repeating: 0, count: n * n)
    LWORK = Int32(WORK[0])
    LIWORK = Int32(IWORK[0])
    WORK = [Float](repeating: 0, count: Int(LWORK))
    IWORK = [Int32](repeating: 0, count: Int(LIWORK))
    sstevd_(
      &JOBZ, &N, &D, &E, &Z, &LDZ, &WORK, &LWORK, &IWORK, &LIWORK, &INFO)
    print("INFO:", INFO)
    print("D:", D)
    print("E:", E)
    print("Z:", Z)
    
    // Check that the eigenvectors produce the eigenvalues.
    let HΨ = Self.matrixMultiply(
      matrixA: originalMatrixT, matrixB: Z, transposeB: true, n: n)
    for electronID in 0..<n {
      var actualE: Float = .zero
      for cellID in 0..<n {
        let address = cellID * n + electronID
        actualE += HΨ[address] * HΨ[address]
      }
      actualE.formSquareRoot()
      
      let expectedE = D[electronID]
      XCTAssertEqual(actualE, expectedE.magnitude, accuracy: 1e-5)
    }
  }
}
