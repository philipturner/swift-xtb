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
  //
  // No intermediate householder reflectors are returned, as this isn't
  // the best algorithm for finding eigenvectors. When debugging, use
  // tridiagonalization as a reference for eigenvalues, and the power method as
  // a reference for eigenvectors.
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
    let nb: Int = 3
    
    var lapackA = [Float](repeating: 0, count: n * n)
    var lapackT = [Float](repeating: 0, count: n * n)
    do {
      for i in originalMatrixA.indices {
        lapackA[i] = originalMatrixA[i]
      }
      lapackA = Self.transpose(matrix: lapackA, n: n)
      
      var M = Int32(n)
      var N = Int32(n)
      var LDA = Int32(n)
      var LDT = Int32(n)
      var INFO = Int32(0)
      sgeqrt2_(&M, &N, &lapackA, &LDA, &lapackT, &LDT, &INFO)
      XCTAssertEqual(INFO, 0, "Received LAPACK error code: \(INFO)")
      
      // Return to row-major format.
      lapackA = Self.transpose(matrix: lapackA, n: n)
      lapackT = Self.transpose(matrix: lapackT, n: n)
      
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
      for rowID in 0..<n {
        for columnID in 0..<n {
          let address = rowID * n + columnID
          let value = lapackT[address]
          print(value, terminator: ", ")
        }
        print()
      }
    }
    
    var currentMatrixA = originalMatrixA
    var currentMatrixT = [Float](repeating: 0, count: n * n)
    var currentMatrixV = [Float](repeating: 0, count: n * n)
    
    var blockStart: Int = 0
    while blockStart < n {
      let blockEnd = min(blockStart + nb, n)
      defer { blockStart += nb }
      
      // Load the panel into the cache. Note that, until some debugging is done,
      // the panel memory allocation will be excessively large.
      var panelMatrixA = [Float](repeating: 0, count: n * n)
      panelMatrixA = currentMatrixA
      
      // Generate the reflectors.
      for transformID in blockStart..<blockEnd {
        // Load the column into the cache.
        var v = [Float](repeating: 0, count: n)
        var columnNorm: Float = .zero
        for rowID in transformID..<n {
          let address = rowID * n + transformID
          let entry = panelMatrixA[address]
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
        
        // Store the column in main memory.
        for rowID in 0..<n {
          let address = rowID * n + transformID
          currentMatrixV[address] = v[rowID]
        }
        
        // Apply H(i) to A(I:M,I+1:N) from the left
        let τ = (newSubdiagonal - oldSubdiagonal) / newSubdiagonal
        for columnID in blockStart..<blockEnd {
          var dotProduct: Float = .zero
          for rowID in 0..<n {
            let matrixAddress = rowID * n + columnID
            let matrixValue = panelMatrixA[matrixAddress]
            dotProduct += v[rowID] * matrixValue
          }
          for rowID in 0..<n {
            let matrixAddress = rowID * n + columnID
            var matrixValue = panelMatrixA[matrixAddress]
            matrixValue -= τ * v[rowID] * dotProduct
            panelMatrixA[matrixAddress] = matrixValue
          }
        }
      }
      
      // Perform a GEMM with non-square matrices.
      var reflectorDotProducts = [Float](repeating: 0, count: nb * nb)
      for lhsID in blockStart..<blockEnd {
        for rhsID in blockStart..<blockEnd {
          var dotProduct: Float = .zero
          for rowID in 0..<n {
            let lhsAddress = rowID * n + lhsID
            let rhsAddress = rowID * n + rhsID
            let lhsValue = currentMatrixV[lhsAddress]
            let rhsValue = currentMatrixV[rhsAddress]
            dotProduct += lhsValue * rhsValue
          }
          let address = (lhsID - blockStart) * nb + (rhsID - blockStart)
          reflectorDotProducts[address] = dotProduct
        }
      }
      
      // Generate the diagonal entries.
      for transformID in blockStart..<blockEnd {
        let diagonalAddress = transformID * n + transformID
        let dotProductAddress = (
          transformID - blockStart) * nb + (transformID - blockStart)
        let dotProduct = reflectorDotProducts[dotProductAddress]
        currentMatrixT[diagonalAddress] = 2 / dotProduct
      }
      
      // Generate the T matrix.
      for transformID in blockStart..<blockEnd {
        // Load the column into the cache.
        var t = [Float](repeating: 0, count: nb)
        for rowID in blockStart..<transformID {
          let address = (rowID - blockStart) * nb + (transformID - blockStart)
          t[rowID - blockStart] = reflectorDotProducts[address]
        }
        
        // This GEMV operation could be transformed into a panel GEMM,
        // similarly to the method for accelerating classical Gram-Schmidt.
        var tt = [Float](repeating: 0, count: nb)
        let τ = currentMatrixT[transformID * n + transformID]
        for rowID in blockStart..<blockEnd {
          // Multiply with the preceding submatrix.
          var dotProduct: Float = .zero
          for columnID in blockStart..<blockEnd {
            let matrixAddress = rowID * n + columnID
            let matrixValue = currentMatrixT[matrixAddress]
            let vectorValue = t[columnID - blockStart]
            dotProduct += matrixValue * vectorValue
          }
          
          // Scale by the value on the diagonal.
          tt[rowID - blockStart] = -τ * dotProduct
        }
        
        // Store the column to main memory.
        for rowID in blockStart..<transformID {
          let address = rowID * n + transformID
          currentMatrixT[address] = tt[rowID - blockStart]
        }
      }
      
      // Update by applying H**T to A(I:M,I+IB:N) from the left
      for columnID in 0..<n {
        // V^H A
        var reflectorContributions = [Float](repeating: 0, count: nb)
        for transformID in blockStart..<blockEnd {
          var dotProduct: Float = .zero
          for rowID in 0..<n {
            let matrixAddress = rowID * n + columnID
            let matrixValue = currentMatrixA[matrixAddress]
            let vectorValue = currentMatrixV[rowID * n + transformID]
            dotProduct += vectorValue * matrixValue
          }
          reflectorContributions[transformID - blockStart] = dotProduct
        }
        
        // T^H (V^H A)
        var TVA = [Float](repeating: 0, count: nb)
        for rowID in blockStart..<blockEnd {
          var dotProduct: Float = .zero
          for columnID in blockStart..<blockEnd {
            let addressT = columnID * n + rowID
            let valueVA = reflectorContributions[columnID - blockStart]
            dotProduct += currentMatrixT[addressT] * valueVA
          }
          TVA[rowID - blockStart] = dotProduct
        }
        
        // V T^H V^H A
        for transformID in blockStart..<blockEnd {
          let dotProduct = TVA[transformID - blockStart]
          for rowID in 0..<n {
            let matrixAddress = rowID * n + columnID
            var matrixValue = currentMatrixA[matrixAddress]
            let vectorValue = currentMatrixV[rowID * n + transformID]
            matrixValue -= dotProduct * vectorValue
            currentMatrixA[matrixAddress] = matrixValue
          }
        }
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
    for rowID in 0..<n {
      for columnID in 0..<n {
        let address = rowID * n + columnID
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
    
    //  LAPACK A
    //  -69.64926, -23.99164, -9.865144, 46.70545, -7.652917, -1.6654878, -9.477488,
    //  0.07827864, -6.594033, -6.417779, -11.78948, -4.3060145, -3.203233, -0.70056605,
    //  0.02609288, -0.15850224, -1.7294693, -4.1016426, -5.185313, -7.5237455, -8.966424,
    //  0.001304644, 0.014304366, 0.16977316, 5.028506, 2.067993, 5.6404796, 6.2869377,
    //  -0.001304644, 0.021262798, -0.3551974, -0.52599597, -3.4273968, -4.896879, -5.8522215,
    //  -0.01304644, -0.14304365, -0.68090385, -0.5504439, -0.5724207, -0.802663, -1.7332793,
    //  0.90020436, -0.35554764, -0.37410447, -0.11177056, 0.5440127, 0.55618095, 0.10665422,
    //
    //  LAPACK T
    //  1.1005036, 0.45809716, 0.34207067, -0.2758081, -0.6424131, -0.6396352, 0.0,
    //  0.0, 1.7055311, -0.12954077, -0.13469173, 0.117625, 0.6343695, 0.0,
    //  0.0, 0.0, 1.1372862, -1.1046314, -0.13303497, 0.45582017, 0.0,
    //  0.0, 0.0, 0.0, 1.2561607, 0.42044175, 1.3487611, 0.0,
    //  0.0, 0.0, 0.0, 0.0, 1.2318188, 0.5077497, 0.0,
    //  0.0, 0.0, 0.0, 0.0, 0.0, 1.5274904, 0.0,
    //  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    //
    //  Custom A
    //  -69.64927, -23.99164, -9.865144, 46.705444, -7.652917, -1.665488, -9.477488,
    //  1.4305115e-06, -6.5940332, -6.4177775, -11.789484, -4.3060136, -3.203233, -0.7005658,
    //  -4.3554346e-06, 0.0, -1.7294698, -4.1016445, -5.185313, -7.5237455, -8.966423,
    //  1.8115559e-06, 0.0, -6.4430594e-08, 5.0285063, 2.0679917, 5.640477, 6.2869344,
    //  5.5420805e-06, 0.0, -6.3509935e-07, -2.3841858e-07, -3.4273992, -4.896884, -5.8522277,
    //  4.6274045e-06, 0.0, -3.99526e-07, -1.0194293e-07, -2.3841858e-07, -0.80266285, -1.7332785,
    //  5.158254e-07, 0.0, -3.620238e-08, 1.2970284e-07, -1.0581476e-07, 4.172325e-07, -0.10665393,
    //
    //  Custom T
    //  1.1005036, 0.45809713, 0.3420709, 0.0, 0.0, 0.0, 0.0,
    //  0.0, 1.705531, -0.12954094, 0.0, 0.0, 0.0, 0.0,
    //  0.0, 0.0, 1.1372862, 0.0, 0.0, 0.0, 0.0,
    //  0.0, 0.0, 0.0, 1.2561606, 0.42044175, 1.348761, 0.0,
    //  0.0, 0.0, 0.0, 0.0, 1.2318188, 0.5077489, 0.0,
    //  0.0, 0.0, 0.0, 0.0, 0.0, 1.5274906, 0.0,
    //  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0,
    //
    //  Custom V
    //  1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    //  0.07827864, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    //  0.02609288, -0.15850224, 1.0, 0.0, 0.0, 0.0, 0.0,
    //  0.001304644, 0.014304366, 0.16977315, 1.0, 0.0, 0.0, 0.0,
    //  -0.001304644, 0.021262798, -0.35519737, -0.52599597, 1.0, 0.0, 0.0,
    //  -0.01304644, -0.14304365, -0.6809038, -0.55044407, -0.57242036, 1.0, 0.0,
    //  0.90020436, -0.35554767, -0.3741047, -0.111770384, 0.54401314, 0.5561807, 1.0,
  }
}
