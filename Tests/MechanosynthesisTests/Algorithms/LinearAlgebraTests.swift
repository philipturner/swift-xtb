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
    let nb: Int = 4
    
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
    
    // TODO: Wrap this into a utility function that accepts an 'n by nb' panel.
    // func panelFactorize(matrix:n:nb:) -> (V: [Float], T: [Float])
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
    
    // Compare to LAPACK results.
    
    // TODO:
    // - Embed the raw numbers generated by LAPACK in assert statements.
    // - Factorize all of the panels of the matrix in Swift.
    // - Remove the dependency on the LAPACK function that requires Swift 5.9.
  }
}

#if false
func referenceCode() {
  // Execute the compact WY transform. This consists of several matrix
  // multiplications, which are split up for clarity.
  print()
  print("Compact WY Transform")
  
  // Operation 1: V T
  var VT = [Float](repeating: 0, count: n * nb)
  for rowID in 0..<n {
    for columnID in blockStart..<blockEnd {
      var dotProduct: Float = .zero
      for innerLoopID in blockStart..<blockEnd {
        let valueV = currentMatrixV[rowID * n + innerLoopID]
        let valueT = currentMatrixT[innerLoopID * n + columnID]
        dotProduct += valueV * valueT
      }
      let addressVT = rowID * nb + (columnID - blockStart)
      VT[addressVT] = dotProduct
    }
  }
  print("VT", VT)
  
  // Operation 2: A V T
  var X = [Float](repeating: 0, count: n * nb)
  for rowID in 0..<n {
    for columnID in blockStart..<blockEnd {
      var dotProduct: Float = .zero
      for innerLoopID in 0..<n {
        let addressA = rowID * n + innerLoopID
        let addressVT = innerLoopID * nb + (columnID - blockStart)
        dotProduct += currentMatrixA[addressA] * VT[addressVT]
      }
      let addressX = rowID * nb + (columnID - blockStart)
      X[addressX] = dotProduct
    }
  }
  print("X", X)
  
  // Operation 3: V^H X
  var VX = [Float](repeating: 0, count: nb * nb)
  for lhsID in blockStart..<blockEnd {
    for rhsID in blockStart..<blockEnd {
      var dotProduct: Float = .zero
      for innerLoopID in 0..<n {
        let addressV = innerLoopID * n + lhsID
        let addressX = innerLoopID * nb + (rhsID - blockStart)
        dotProduct += currentMatrixV[addressV] * X[addressX]
      }
      let addressVX = (lhsID - blockStart) * nb + (rhsID - blockStart)
      VX[addressVX] = dotProduct
    }
  }
  print("VX", VX)
  
  // Operation 4: T^H V^H X
  var TVX = [Float](repeating: 0, count: nb * nb)
  for lhsID in blockStart..<blockEnd {
    for rhsID in blockStart..<blockEnd {
      var dotProduct: Float = .zero
      for innerLoopID in blockStart..<blockEnd {
        let addressT = innerLoopID * n + lhsID
        let addressVX = (
          innerLoopID - blockStart) * nb + (rhsID - blockStart)
        dotProduct += currentMatrixT[addressT] * VX[addressVX]
      }
      let addressTVX = (lhsID - blockStart) * nb + (rhsID - blockStart)
      TVX[addressTVX] = dotProduct
    }
  }
  print("TVX", TVX)
  
  // Operation 5: X - (1 / 2) V (T^H V^H X)
  var W = [Float](repeating: 0, count: n * nb)
  for rowID in 0..<n {
    for columnID in blockStart..<blockEnd {
      var dotProduct: Float = .zero
      for innerLoopID in blockStart..<blockEnd {
        let addressV = rowID * n + innerLoopID
        let addressTVX = (
          innerLoopID - blockStart) * nb + (columnID - blockStart)
        dotProduct += currentMatrixV[addressV] * TVX[addressTVX]
      }
      let address = rowID * nb + (columnID - blockStart)
      W[address] = X[address] - 0.5 * dotProduct
    }
  }
  print("W", W)
  
  // Operation 6: A - W V^H
  for rowID in 0..<n {
    for columnID in 0..<n {
      var dotProduct: Float = .zero
      for innerLoopID in blockStart..<blockEnd {
        let addressW = rowID * nb + (innerLoopID - blockStart)
        let addressV = columnID * n + innerLoopID
        dotProduct += W[addressW] * currentMatrixV[addressV]
      }
      let address = rowID * n + columnID
      let matrixValue = currentMatrixA[address]
      currentMatrixA[address] = matrixValue - dotProduct
    }
  }
  print("A", currentMatrixA)
  
  // Copy the panel into the original matrix.
  for rowID in 0..<n {
    for columnID in blockStart..<blockEnd {
      let address = rowID * n + columnID
      currentMatrixA[address] = panelMatrixA[address]
    }
  }
  
  // TODO: Omit when the WY transform is finished.
  //      currentMatrixA = panelMatrixA
}
#endif
