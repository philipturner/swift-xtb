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
  
  // Diagonalizes a tridiagonal matrix with LAPACK divide-and-conquer.
  // Returns the eigenvectors in column-major format, as output by LAPACK.
  static func divideAndConquer(matrix: [Float], n: Int) -> (
    eigenvalues: [Float], eigenvectors: [Float]
  ) {
    // Store the tridiagonal matrix in a compact form.
    var D = [Float](repeating: 0, count: n)
    var E = [Float](repeating: 0, count: n - 1)
    for diagonalID in 0..<n {
      let matrixAddress = diagonalID * n + diagonalID
      let vectorAddress = diagonalID
      D[vectorAddress] = matrix[matrixAddress]
    }
    for subDiagonalID in 0..<n - 1 {
      let rowID = subDiagonalID
      let columnID = subDiagonalID + 1
      let matrixAddress = rowID * n + columnID
      let vectorAddress = rowID
      E[vectorAddress] = matrix[matrixAddress]
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
    
    // Return the eigenpairs.
    return (eigenvalues: D, eigenvectors: Z)
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
      let (D, Z) = Self.divideAndConquer(matrix: originalMatrixT, n: n)
      
      // Check that the eigenvectors produce the eigenvalues.
      let HΨ = Self.matrixMultiply(
        matrixA: originalMatrixT, 
        matrixB: Z,
        transposeB: true, n: n)
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
  
  // Test the two-stage process for tridiagonalizing a matrix, and the
  // aggregation of the Householder reflectors.
  func testTwoStageTridiagonalization() throws {
    // TODO:
    // - Find an algorithm for programming the order of bulge chases.
    // - Test the full-stack diagonalizer against adversarial edge cases.
    // - Run benchmarks, starting with the unoptimized single-core CPU code.
    //   This marks the beginning of the optimization process, which should
    //   end with GPU acceleration.
    var originalMatrixA: [Float] = [
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
    
    // Make the matrix symmetric.
    originalMatrixA = Self.matrixMultiply(
      matrixA: originalMatrixA,
      matrixB: originalMatrixA,
      transposeB: true, n: n)
    
    // Allocate main memory allocations.
    var currentMatrixA = originalMatrixA
    var currentReflectors = [Float](repeating: 0, count: n * n)
    
    // Reduce the matrix to band form, and collect up the reflectors.
    var blockStart: Int = 0
    while blockStart < n - nb {
      // Adjust the loop end, to account for the factorization band offset.
      let blockEnd = min(blockStart + nb, n - nb)
      defer { blockStart += nb }
      
      // Load to panel into the cache, isolating mutations to the matrix A.
      var panel = [Float](repeating: 0, count: nb * n)
      for rowID in blockStart..<blockEnd {
        for columnID in 0..<n {
          let matrixAddress = rowID * n + columnID
          let panelAddress = (rowID - blockStart) * n + columnID
          panel[panelAddress] = currentMatrixA[matrixAddress]
        }
      }
      
      // Allocate cache memory for the reflectors.
      var panelReflectors = [Float](repeating: 0, count: nb * n)
      
      // Generate the reflectors.
      for reflectorID in blockStart..<blockEnd {
        // Factor starting at an offset from the diagonal.
        let bandOffset = reflectorID + nb
        
        // Load the row into the cache.
        var vector = [Float](repeating: 0, count: n)
        for elementID in 0..<n {
          let address = (reflectorID - blockStart) * n + elementID
          vector[elementID] = panel[address]
        }
        
        // Apply preceding reflectors (from this panel) to the column.
        for previousReflectorID in blockStart..<reflectorID {
          // Load the reflector into the cache.
          var reflector = [Float](repeating: 0, count: n)
          for elementID in 0..<n {
            let address = (previousReflectorID - blockStart) * n + elementID
            reflector[elementID] = panelReflectors[address]
          }
          
          // Apply the reflector.
          var dotProduct: Float = .zero
          for elementID in 0..<n {
            dotProduct += reflector[elementID] * vector[elementID]
          }
          for elementID in 0..<n {
            vector[elementID] -= reflector[elementID] * dotProduct
          }
        }
        
        // Zero out the elements above the band offset.
        for elementID in 0..<bandOffset {
          vector[elementID] = 0
        }
        
        // Take the norm of the vector.
        var norm: Float = .zero
        for elementID in 0..<n {
          norm += vector[elementID] * vector[elementID]
        }
        norm.formSquareRoot()
        
        // Modify the vector, turning it into a reflector.
        let oldSubdiagonal = vector[bandOffset]
        let newSubdiagonal = norm * Float((oldSubdiagonal >= 0) ? -1 : 1)
        let tau = (newSubdiagonal - oldSubdiagonal) / newSubdiagonal
        for elementID in 0..<n {
          var element = vector[elementID]
          if elementID == bandOffset {
            element = 1
          } else {
            element /= oldSubdiagonal - newSubdiagonal
          }
          element *= tau.squareRoot()
          vector[elementID] = element
        }
        
        // Store the reflector to the cache.
        for elementID in 0..<n {
          let address = (reflectorID - blockStart) * n + elementID
          panelReflectors[address] = vector[elementID]
        }
      }
      
      // Apply the reflectors to the matrix, from both sides.
      for directionID in 0..<2 {
        for vectorID in 0..<n {
          var vector = [Float](repeating: 0, count: n)
          if directionID == 0 {
            // Load the row into the cache.
            for elementID in 0..<n {
              let address = vectorID * n + elementID
              vector[elementID] = currentMatrixA[address]
            }
          } else {
            // Load the column into the cache.
            for elementID in 0..<n {
              let address = elementID * n + vectorID
              vector[elementID] = currentMatrixA[address]
            }
          }
          
          for reflectorID in blockStart..<blockEnd {
            // Load the reflector into the cache.
            var reflector = [Float](repeating: 0, count: n)
            for elementID in 0..<n {
              let address = (reflectorID - blockStart) * n + elementID
              reflector[elementID] = panelReflectors[address]
            }
            
            // Apply the reflector.
            var dotProduct: Float = .zero
            for elementID in 0..<n {
              dotProduct += reflector[elementID] * vector[elementID]
            }
            for elementID in 0..<n {
              vector[elementID] -= reflector[elementID] * dotProduct
            }
          }
          
          if directionID == 0 {
            // Store the row to main memory.
            for elementID in 0..<n {
              let address = vectorID * n + elementID
              currentMatrixA[address] = vector[elementID]
            }
          } else {
            // Store the column to main memory.
            for elementID in 0..<n {
              let address = elementID * n + vectorID
              currentMatrixA[address] = vector[elementID]
            }
          }
        }
      }
      
      // Store the reflectors to main memory.
      for reflectorID in blockStart..<blockEnd {
        for elementID in 0..<n {
          let cacheAddress = (reflectorID - blockStart) * n + elementID
          let memoryAddress = reflectorID * n + elementID
          currentReflectors[memoryAddress] = panelReflectors[cacheAddress]
        }
      }
    }
    
    print()
    print("Matrix A")
    for rowID in 0..<n {
      for columnID in 0..<n {
        let address = rowID * n + columnID
        var value = currentMatrixA[address]
        if value.magnitude < 1e-3 {
          value = 0
        }
        print(value, terminator: ", ")
      }
      print()
    }
    
    // MARK: - Bulge Chasing
    
    var bulgeReflectors: [[Float]] = []
    
    #if false
    // nb = 3
    chaseBulge(reflectorID: 0, diagonalOffset: 2)
    chaseBulge(reflectorID: 2, diagonalOffset: 3)
    chaseBulge(reflectorID: 1, diagonalOffset: 2)
    chaseBulge(reflectorID: 2, diagonalOffset: 2)
    chaseBulge(reflectorID: 3, diagonalOffset: 2)
    
    chaseBulge(reflectorID: 0, diagonalOffset: 1)
    chaseBulge(reflectorID: 1, diagonalOffset: 2)
    chaseBulge(reflectorID: 3, diagonalOffset: 2)
    
    chaseBulge(reflectorID: 1, diagonalOffset: 1)
    chaseBulge(reflectorID: 2, diagonalOffset: 2)
    
    chaseBulge(reflectorID: 2, diagonalOffset: 1)
    chaseBulge(reflectorID: 3, diagonalOffset: 2)
    
    chaseBulge(reflectorID: 3, diagonalOffset: 1)
    chaseBulge(reflectorID: 4, diagonalOffset: 1)
    #endif
    
    #if false
    // nb = 4
    chaseBulge(reflectorID: 0, diagonalOffset: 3)
    chaseBulge(reflectorID: 1, diagonalOffset: 3)
    chaseBulge(reflectorID: 2, diagonalOffset: 3)
    
    chaseBulge(reflectorID: 0, diagonalOffset: 2)
    chaseBulge(reflectorID: 2, diagonalOffset: 3)
    chaseBulge(reflectorID: 1, diagonalOffset: 2)
    chaseBulge(reflectorID: 2, diagonalOffset: 2)
    chaseBulge(reflectorID: 3, diagonalOffset: 2)
    
    chaseBulge(reflectorID: 0, diagonalOffset: 1)
    chaseBulge(reflectorID: 1, diagonalOffset: 2)
    chaseBulge(reflectorID: 3, diagonalOffset: 2)
    
    chaseBulge(reflectorID: 1, diagonalOffset: 1)
    chaseBulge(reflectorID: 2, diagonalOffset: 2)
    
    chaseBulge(reflectorID: 2, diagonalOffset: 1)
    chaseBulge(reflectorID: 3, diagonalOffset: 2)
    
    chaseBulge(reflectorID: 3, diagonalOffset: 1)
    chaseBulge(reflectorID: 4, diagonalOffset: 1)
    #endif
    
    var bulgeMask = [Int](repeating: 0, count: n * n)
    for rowID in 0..<n {
      bulgeMask[rowID * n + rowID] = 1
      
      for subDiagonalID in 1...nb {
        guard rowID + subDiagonalID < n else {
          continue
        }
        bulgeMask[(rowID + subDiagonalID) * n + rowID] = 1
      }
    }
    
    print()
    print("bulge mask")
    for rowID in 0..<n {
      var maskRow = [Int](repeating: 0, count: n)
      for columnID in 0..<n {
        let address = rowID * n + columnID
        maskRow[columnID] = bulgeMask[address]
      }
      print(maskRow)
    }
    
    updateBulgeMask(reflectorID: 0, diagonalOffset: 2)
    chaseBulge(reflectorID: 0, diagonalOffset: 2)
    
    updateBulgeMask(reflectorID: 2, diagonalOffset: 3)
    chaseBulge(reflectorID: 2, diagonalOffset: 3)
    
    // TODO: An algorithm that searches for the diagonal ID of the next bulge.
    // Then, an algorithm to chase the bulge down, until that diagonal no
    // longer contains any bulges.
    
    func updateBulgeMask(reflectorID: Int, diagonalOffset: Int) {
      for subDiagonalID in min(diagonalOffset + 1, n)...n {
        let rowID = reflectorID + subDiagonalID
        var address = rowID * n + reflectorID
        if address < n * n {
          bulgeMask[address] = 0
        }
      }
      
      // We're assuming the Householder reflectors are effectively Givens
      // rotations.
      for rowID in 0..<n {
        var containsOne = false
        for subDiagonalID in diagonalOffset...(diagonalOffset + 1) {
          let columnID = reflectorID + subDiagonalID
          var address = rowID * n + columnID
          if address < n * n {
            if bulgeMask[address] > 0 {
              containsOne = true
            }
          }
        }
        if containsOne {
          for subDiagonalID in diagonalOffset...(diagonalOffset + 1) {
            let columnID = reflectorID + subDiagonalID
            var address = rowID * n + columnID
            if address < n * n {
              bulgeMask[address] = 2
            }
          }
        }
      }
      
      print()
      print("bulge mask")
      for rowID in 0..<n {
        var maskRow = [Int](repeating: 0, count: n)
        for columnID in 0..<n {
          let address = rowID * n + columnID
          maskRow[columnID] = bulgeMask[address]
        }
        print(maskRow)
      }
    }
    
    func chaseBulge(reflectorID: Int, diagonalOffset: Int) {
      let bandOffset = reflectorID + diagonalOffset
      
      // Load the row into the cache.
      var vector = [Float](repeating: 0, count: n)
      for elementID in bandOffset..<n {
        let address = reflectorID * n + elementID
        vector[elementID] = currentMatrixA[address]
      }
      print()
      print("bulge row:", vector)
      
      // Take the norm of the vector.
      var norm: Float = .zero
      for elementID in 0..<n {
        norm += vector[elementID] * vector[elementID]
      }
      norm.formSquareRoot()
      
      // Modify the vector, turning it into a reflector.
      let oldSubdiagonal = vector[bandOffset]
      let newSubdiagonal = norm * Float((oldSubdiagonal >= 0) ? -1 : 1)
      let tau = (newSubdiagonal - oldSubdiagonal) / newSubdiagonal
      for elementID in 0..<n {
        var element = vector[elementID]
        if elementID == bandOffset {
          element = 1
        } else {
          element /= oldSubdiagonal - newSubdiagonal
        }
        element *= tau.squareRoot()
        vector[elementID] = element
      }
      print("bulge reflector:", vector)
      
      // Apply the reflector to the matrix, from both sides.
      let reflector = vector
      for directionID in 0..<2 {
        for vectorID in 0..<n {
          var vector = [Float](repeating: 0, count: n)
          if directionID == 0 {
            // Load the row into the cache.
            for elementID in 0..<n {
              let address = vectorID * n + elementID
              vector[elementID] = currentMatrixA[address]
            }
          } else {
            // Load the column into the cache.
            for elementID in 0..<n {
              let address = elementID * n + vectorID
              vector[elementID] = currentMatrixA[address]
            }
          }
          
          // Apply the reflector.
          var dotProduct: Float = .zero
          for elementID in 0..<n {
            dotProduct += reflector[elementID] * vector[elementID]
          }
          for elementID in 0..<n {
            vector[elementID] -= reflector[elementID] * dotProduct
          }
          
          if directionID == 0 {
            // Store the row to main memory.
            for elementID in 0..<n {
              let address = vectorID * n + elementID
              currentMatrixA[address] = vector[elementID]
            }
          } else {
            // Store the column to main memory.
            for elementID in 0..<n {
              let address = elementID * n + vectorID
              currentMatrixA[address] = vector[elementID]
            }
          }
        }
      }
      
      for rowID in 0..<n {
        for columnID in 0..<n {
          let address = rowID * n + columnID
          var value = currentMatrixA[address]
          if value.magnitude < 1e-3 {
            value = 0
          }
          print(value, terminator: ", ")
        }
        print()
      }
      
      // Store the reflector to main memory.
      bulgeReflectors.append(reflector)
    }
    
    // MARK: - Validation Testing
    
    return
    
    // Test: Diagonalize the banded matrix with standard techniques. Acquire
    // the eigenvectors, then back-transform them using the reflectors. Ensure
    // they return the same eigenvalue as expected.
    var (eigenvalues, eigenvectors) = Self
      .divideAndConquer(matrix: currentMatrixA, n: n)
    eigenvectors = Self.transpose(matrix: eigenvectors, n: n)
    
    print()
    print("eigenvalues from D&C")
    print(eigenvalues)
    
    // Display the eigenvectors before the transformation.
    print()
    print("eigenvectors before transformation")
    for vectorID in 0..<n {
      var vector = [Float](repeating: 0, count: n)
      var eigenvalue: Float = .zero
      let matrixRow = Int.random(in: 0..<n)
      
      for elementID in 0..<n {
        let vectorAddress = elementID * n + vectorID
        let vectorValue = eigenvectors[vectorAddress]
        vector[elementID] = vectorValue
        
        // Read data from the current storage for matrix A.
        let matrixAddress = matrixRow * n + elementID
        let matrixValue = currentMatrixA[matrixAddress]
        eigenvalue += matrixValue * vectorValue
        
      }
      eigenvalue /= vector[matrixRow]
      print("Ψ[\(eigenvalue)]:", vector)
    }
    
    eigenvectors = Self.transpose(matrix: eigenvectors, n: n)
    
    // Back-transform the eigenvectors.
    print()
    print("back-transforming eigenvectors")
    for vectorID in 0..<n {
      // Load the vector into the cache.
      var vector = [Float](repeating: 0, count: n)
      for elementID in 0..<n {
        let address = vectorID * n + elementID
        vector[elementID] = eigenvectors[address]
      }
      
      for reflectorID in bulgeReflectors.indices.reversed() {
        // Load the reflector into the cache.
        let reflector = bulgeReflectors[reflectorID]
        
        // Apply the reflector.
        var dotProduct: Float = .zero
        for elementID in 0..<n {
          dotProduct += reflector[elementID] * vector[elementID]
        }
        for elementID in 0..<n {
          vector[elementID] -= reflector[elementID] * dotProduct
        }
      }
      
      for reflectorID in (0..<n).reversed() {
        // Load the reflector into the cache.
        var reflector = [Float](repeating: 0, count: n)
        for elementID in 0..<n {
          let address = reflectorID * n + elementID
          reflector[elementID] = currentReflectors[address]
        }
        
        // Apply the reflector.
        var dotProduct: Float = .zero
        for elementID in 0..<n {
          dotProduct += reflector[elementID] * vector[elementID]
        }
        for elementID in 0..<n {
          vector[elementID] -= reflector[elementID] * dotProduct
        }
      }
      print("v[\(vectorID)]", vector)
      
      // Store the vector to main memory.
      for elementID in 0..<n {
        let address = vectorID * n + elementID
        eigenvectors[address] = vector[elementID]
      }
    }
    
    // Display the eigenvectors after the transformation.
    print()
    print("eigenvectors after transformation")
    for vectorID in 0..<n {
      var vector = [Float](repeating: 0, count: n)
      var eigenvalue: Float = .zero
      let matrixRow = Int.random(in: 0..<n)
      for elementID in 0..<n {
        let vectorAddress = vectorID * n + elementID
        let vectorValue = eigenvectors[vectorAddress]
        vector[elementID] = vectorValue
        
        // Read data from the current storage for matrix A.
        let matrixAddress = matrixRow * n + elementID
        let matrixValue = originalMatrixA[matrixAddress]
        eigenvalue += matrixValue * vectorValue
        
      }
      eigenvalue /= vector[matrixRow]
      print("Ψ[\(eigenvalue)]:", vector)
    }
  }
}
