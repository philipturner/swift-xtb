import XCTest
import Accelerate
import Numerics

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
        var dotProduct: Float = .zero
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
          dotProduct += value1 * value2
        }
        matrixC[rowID * n + columnID] = dotProduct
      }
    }
    return matrixC
  }
  
  // Forms an orthogonal basis from a square matrix's columns.
  static func modifiedGramSchmidt(
    matrix originalMatrix: [Float], n: Int
  ) -> [Float] {
    // Operate on the output matrix in-place.
    var matrix = originalMatrix
    
    func normalize(electronID: Int) {
      var norm: Float = .zero
      for cellID in 0..<n {
        let value = matrix[cellID * n + electronID]
        norm += value * value
      }
      
      let normalizationFactor = 1 / norm.squareRoot()
      for cellID in 0..<n {
        var value = matrix[cellID * n + electronID]
        value *= normalizationFactor
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
        var dotProduct: Float = .zero
        for cellID in 0..<n {
          let value1 = matrix[cellID * n + electronID]
          let value2 = matrix[cellID * n + neighborID]
          dotProduct += value1 * value2
        }
        
        // Subtract the parallel component.
        for cellID in 0..<n {
          var value1 = matrix[cellID * n + electronID]
          let value2 = matrix[cellID * n + neighborID]
          value1 -= dotProduct * value2
          matrix[cellID * n + electronID] = value1
        }
      }
      
      // Rescale the orthogonal component to unit vector length.
      normalize(electronID: electronID)
    }
    
    return matrix
  }
  
  // Reduce the matrix to tridiagonal form.
  static func tridiagonalize(
    matrix originalMatrix: [Float],
    n: Int
  ) -> [Float] {
    var currentMatrixA = originalMatrix
    
    // This requires that n > 1.
    for transformID in 0..<n - 2 {
      // Load the column into the cache.
      var V = [Float](repeating: 0, count: n)
      var columnNorm: Float = .zero
      for rowID in (transformID + 1)..<n {
        let address = rowID * n + transformID
        let entry = currentMatrixA[address]
        V[rowID] = entry
        columnNorm += entry * entry
      }
      columnNorm.formSquareRoot()
      
      // Form the 'v' output of Householder(j,x).
      let oldSubdiagonal = V[transformID + 1]
      let newSubdiagonal = columnNorm * Float((oldSubdiagonal >= 0) ? -1 : 1)
      V[transformID + 1] = 1
      for rowID in (transformID + 2)..<n {
        V[rowID] /= oldSubdiagonal - newSubdiagonal
      }
      
      // Form the 'τ' output of Householder(j,x).
      let T = (newSubdiagonal - oldSubdiagonal) / newSubdiagonal
      
      // Operation 1: VT
      var VT = [Float](repeating: 0, count: n)
      for rowID in 0..<n {
        VT[rowID] = V[rowID] * T
      }
      
      // Operation 2: AVT
      var X = [Float](repeating: 0, count: n)
      for rowID in 0..<n {
        var dotProduct: Float = .zero
        for columnID in 0..<n {
          let address = rowID * n + columnID
          dotProduct += currentMatrixA[address] * VT[columnID]
        }
        X[rowID] = dotProduct
      }
      
      // Operation 3: V^H X
      var VX: Float = .zero
      for rowID in 0..<n {
        VX += V[rowID] * X[rowID]
      }
      
      // Operation 4: X - (1 / 2) VT^H (V^H X)
      var W = [Float](repeating: 0, count: n)
      for rowID in 0..<n {
        W[rowID] = X[rowID] - 0.5 * V[rowID] * T * VX
      }
      
      // Operation 5: A - WV^H - VW^H
      for rowID in 0..<n {
        for columnID in 0..<n {
          let address = rowID * n + columnID
          var entry = currentMatrixA[address]
          entry -= W[rowID] * V[columnID]
          entry -= V[rowID] * W[columnID]
          currentMatrixA[address] = entry
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
    // The example sourced from Wikipedia.
    testMatrix([
      4, 1, -2, 2,
      1, 2, 0, 1,
      -2, 0, 3, -2,
      2, 1, -2, -1
    ], n: 4)
    
    var eigenvectors: [Float] = [
      7, 6, 5, 4, 3, 2, 1,
      6, 7, 5, 4, 3, 2, 1,
      2, -1, -1, 1, 2, 6, 8,
      0.1, 0.2, 0.5, 0.5, 0.1, 0.2, 0.5,
      -0.1, 0.2, -0.5, 0.5, -0.1, 0.2, -0.5,
      -1, -2, -3, -5, -7, -9, -10,
      69, 23, 9, -48, 7, 1, 9,
    ]
    eigenvectors = Self.transpose(matrix: eigenvectors, n: 7)
    eigenvectors = Self.modifiedGramSchmidt(matrix: eigenvectors, n: 7)
    
    // Well-conditioned eigenspectra without degenerate clusters.
    testEigenvalues([
      4, 3, 2, 1, 0.1, -1.3, -2.3
    ])
    testEigenvalues([
      4, 3, 2, 1, 0, -1, -2
    ])
    
    // Problematic for the iterative diagonalization experiment.
    testEigenvalues([
      4, 3.01, 3.00, 2.99, 0, -2, -2
    ])
    testEigenvalues([
      4, 3.5, 3, 2, 0, -1, -2
    ])
    testEigenvalues([
      2, 2, 0, -2.99, -3.00, -3.01, -4
    ])
    testEigenvalues([
      2, 1, 0, -2.99, -3.00, -3.01, -4
    ])
    
    // Not problematic for the iterative diagonalization experiment.
    testEigenvalues([
      4, 3.01, 3.00, 2.99, 0, -1, -2
    ])
    testEigenvalues([
      4, 3.5, 3, 2, 0, -2, -2
    ])
    testEigenvalues([
      2, 2, 0, -2, -3, -3.5, -4
    ])
    
    func testEigenvalues(_ eigenvalues: [Float]) {
      var Λ = [Float](repeating: 0, count: 7 * 7)
      for i in 0..<7 {
        let address = i * 7 + i
        let value = eigenvalues[i]
        Λ[address] = value
      }
      let ΣT = eigenvectors
      
      let ΛΣT = Self.matrixMultiply(
        matrixA: Λ, transposeA: false,
        matrixB: ΣT, transposeB: false, n: 7)
      let A = Self.matrixMultiply(
        matrixA: ΣT, transposeA: true,
        matrixB: ΛΣT, transposeB: false, n: 7)
      let AΣT = Self.matrixMultiply(
        matrixA: A, transposeA: false,
        matrixB: ΣT, transposeB: true, n: 7)
      
      for electronID in 0..<7 {
        let expectedEigenvalue = eigenvalues[electronID]
        var actualEigenvalue: Float = .zero
        for cellID in 0..<7 {
          let value = AΣT[cellID * 7 + electronID]
          actualEigenvalue += value * value
        }
        actualEigenvalue.formSquareRoot()
        XCTAssertLessThan(
          expectedEigenvalue.magnitude - actualEigenvalue.magnitude, 1e-4)
      }
      
      testMatrix(A, n: 7)
    }
    
    func testMatrix(_ originalMatrixA: [Float], n: Int) {
      let originalMatrixT = Self.tridiagonalize(matrix: originalMatrixA, n: n)
      
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
      
      // Call into LAPACK.
      var Z = [Float](repeating: 0, count: n * n)
      LWORK = Int32(WORK[0])
      LIWORK = Int32(IWORK[0])
      WORK = [Float](repeating: 0, count: Int(LWORK))
      IWORK = [Int32](repeating: 0, count: Int(LIWORK))
      sstevd_(
        &JOBZ, &N, &D, &E, &Z, &LDZ, &WORK, &LWORK, &IWORK, &LIWORK, &INFO)
      
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
  
  // Reproduce QR panel factorization (xGEQRT) from the Dongarra 2013 paper.
  func testPanelFactorization() {
    let originalMatrixA: [Float] = [
      7, 6, 5, 4, 3, 2, 1,
      6, 7, 5, 4, 3, 2, 1,
      2, -1, -1, 1, 2, 6, 8,
      0.1, 0.2, 0.5, 0.5, 0.1, 0.2, 0.5,
      -0.1, 0.2, -0.5, 0.5, -0.1, 0.2, -0.5,
      -1, -2, -3, -5, -7, -9, -10,
      69, 23, 9, -48, 7, 1, 9,
    ]
    let n: Int = 7
    let blockSize: Int = 3
    
    var lapackA = [Float](repeating: 0, count: n * n)
    var lapackT = [Float](repeating: 0, count: blockSize * blockSize)
    do {
      for i in originalMatrixA.indices {
        lapackA[i] = originalMatrixA[i]
      }
      lapackA = Self.transpose(matrix: lapackA, n: n)
      
      var M = Int32(n)
      var N = Int32(blockSize)
      var LDA = Int32(n)
      var LDT = Int32(blockSize)
      var INFO = Int32(0)
      sgeqrt2_(&M, &N, &lapackA, &LDA, &lapackT, &LDT, &INFO)
      XCTAssertEqual(INFO, 0, "Received LAPACK error code: \(INFO)")
      
      // Return to row-major format.
      lapackA = Self.transpose(matrix: lapackA, n: n)
      lapackT = Self.transpose(matrix: lapackT, n: blockSize)
      
      print()
      print("LAPACK A")
      for rowID in 0..<n {
        for columnID in 0..<n {
          let address = rowID * n + columnID
          let value = lapackA[address]
          print(value, terminator: ", ")
        }
        print()
      }
      
      print()
      print("LAPACK T")
      for rowID in 0..<blockSize {
        for columnID in 0..<blockSize {
          let address = rowID * blockSize + columnID
          let value = lapackT[address]
          print(value, terminator: ", ")
        }
        print()
      }
    }
    
    
    
    var currentMatrixA = originalMatrixA
    var currentMatrixT = [Float](repeating: 0, count: blockSize * blockSize)
    var currentMatrixV = [Float](repeating: 0, count: n * n)
    
    // Generate the reflectors.
    for transformID in 0..<blockSize {
      // Load the column into the cache.
      var v = [Float](repeating: 0, count: n)
      var columnNorm: Float = .zero
      for rowID in transformID..<n {
        let address = rowID * n + transformID
        let entry = currentMatrixA[address]
        v[rowID] = entry
        columnNorm += entry * entry
      }
      columnNorm.formSquareRoot()
      
      // Generate elem. refl. H(i) to annihilate A(i+1:m,i), tau(I) -> T(I,1)
      let oldSubdiagonal = v[transformID]
      let newSubdiagonal = columnNorm * Float((oldSubdiagonal >= 0) ? -1 : 1)
      v[transformID] = 1
      for rowID in (transformID + 1)..<n {
        v[rowID] /= oldSubdiagonal - newSubdiagonal
      }
      let τ = (newSubdiagonal - oldSubdiagonal) / newSubdiagonal
      
      // Store the column in main memory.
      for rowID in 0..<n {
        let address = rowID * n + transformID
        currentMatrixV[address] = v[rowID]
      }
      currentMatrixT[transformID * blockSize + transformID] = τ
      
      // Apply H(i) to A(I:M,I+1:N) from the left
      for columnID in 0..<blockSize {
        var dotProduct: Float = .zero
        for rowID in 0..<n {
          let matrixAddress = rowID * n + columnID
          let matrixValue = currentMatrixA[matrixAddress]
          dotProduct += v[rowID] * matrixValue
        }
        for rowID in 0..<n {
          let matrixAddress = rowID * n + columnID
          var matrixValue = currentMatrixA[matrixAddress]
          matrixValue -= τ * v[rowID] * dotProduct
          currentMatrixA[matrixAddress] = matrixValue
        }
      }
    }
    
    // Perform a GEMM with non-square matrices.
    var reflectorDotProducts = [Float](
      repeating: 0, count: blockSize * blockSize)
    for lhsID in 0..<blockSize {
      for rhsID in 0..<blockSize {
        var dotProduct: Float = .zero
        for rowID in 0..<n {
          let lhsAddress = rowID * n + lhsID
          let rhsAddress = rowID * n + rhsID
          let lhsValue = currentMatrixV[lhsAddress]
          let rhsValue = currentMatrixV[rhsAddress]
          dotProduct += lhsValue * rhsValue
        }
        let address = lhsID * blockSize + rhsID
        reflectorDotProducts[address] = dotProduct
      }
    }
    
    // Generate the T matrix.
    for transformID in 1..<blockSize {
      // Load the column into the cache.
      var t = [Float](repeating: 0, count: blockSize)
      for rowID in 0..<transformID {
        let address = rowID * blockSize + transformID
        t[rowID] = reflectorDotProducts[address]
      }
      
      // This GEMV operation could be transformed into a panel GEMM,
      // similar to the method used for classical Gram-Schmidt.
      var tt = [Float](repeating: 0, count: blockSize)
      let tau = currentMatrixT[transformID * blockSize + transformID]
      for rowID in 0..<blockSize {
        // Multiply with the preceding submatrix.
        var dotProduct: Float = .zero
        for columnID in 0..<blockSize {
          let matrixAddress = rowID * blockSize + columnID
          let matrixValue = currentMatrixT[matrixAddress]
          let vectorValue = t[columnID]
          dotProduct += matrixValue * vectorValue
        }
        
        // Scale by the value on the diagonal.
        tt[rowID] = -tau * dotProduct
      }
      
      // Store the column to main memory.
      for rowID in 0..<transformID {
        let address = rowID * blockSize + transformID
        currentMatrixT[address] = tt[rowID]
      }
    }
    
    print()
    print("Custom A")
    for rowID in 0..<n {
      for columnID in 0..<n {
        let address = rowID * n + columnID
        let value = currentMatrixA[address]
        print(value, terminator: ", ")
      }
      print()
    }
    
    print()
    print("Custom T")
    for rowID in 0..<blockSize {
      for columnID in 0..<blockSize {
        let address = rowID * blockSize + columnID
        let value = currentMatrixT[address]
        print(value, terminator: ", ")
      }
      print()
    }
    
    print()
    print("Custom V")
    for rowID in 0..<n {
      for columnID in 0..<n {
        let address = rowID * n + columnID
        let value = currentMatrixV[address]
        print(value, terminator: ", ")
      }
      print()
    }
    
    // Compare to LAPACK results.
    
    // TODO:
    // - Embed the raw numbers generated by LAPACK in assert statements.
    // - Factorize all of the panels of the matrix.
    // - Remove the dependency on the LAPACK function that requires Swift 5.9.
  }
}
