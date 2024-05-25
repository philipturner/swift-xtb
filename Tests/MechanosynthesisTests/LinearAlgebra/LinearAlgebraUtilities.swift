//
//  LinearAlgebraUtilities.swift
//
//
//  Created by Philip Turner on 5/8/24.
//

import Accelerate

struct LinearAlgebraUtilities {
  // Multiplies two square matrices.
  static func matrixMultiply(
    matrixA: [Float], transposeA: Bool = false,
    matrixB: [Float], transposeB: Bool = false,
    n: Int
  ) -> [Float] {
    var matrixC = [Float](repeating: 0, count: n * n)
    
#if false
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
#else
    var TRANSA: CChar
    var TRANSB: CChar
    if transposeB {
      TRANSA = CChar(Character("T").asciiValue!)
    } else {
      TRANSA = CChar(Character("N").asciiValue!)
    }
    if transposeA {
      TRANSB = CChar(Character("T").asciiValue!)
    } else {
      TRANSB = CChar(Character("N").asciiValue!)
    }
    
    var M: Int32 = Int32(n)
    var N: Int32 = Int32(n)
    var K: Int32 = Int32(n)
    var ALPHA: Float = 1
    var LDA: Int32 = Int32(n)
    var BETA: Float = 0
    var LDB: Int32 = Int32(n)
    var LDC: Int32 = Int32(n)
    matrixA.withContiguousStorageIfAvailable {
      let B = UnsafeMutablePointer(mutating: $0.baseAddress!)
      matrixB.withContiguousStorageIfAvailable {
        let A = UnsafeMutablePointer(mutating: $0.baseAddress!)
        matrixC.withContiguousMutableStorageIfAvailable {
          let C = $0.baseAddress!
          sgemm_(
            &TRANSA, // TRANSA
            &TRANSB, // TRANSB
            &M, // M
            &N, // N
            &K, // K
            &ALPHA, // ALPHA
            A, // A
            &LDA, // LDA
            B, // B
            &LDB, // LDB
            &BETA, // BETA
            C, // C
            &LDC // LDC
          )
        }
      }
    }
#endif
    
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
  
  // Reference implementation for diagonalization in LAPACK.
  static func diagonalize(matrix: [Float], n: Int) -> (
    eigenvalues: [Float], eigenvectors: [Float]
  ) {
    var JOBZ = CChar(Character("V").asciiValue!)
    var UPLO = CChar(Character("L").asciiValue!)
    var N: Int32 = Int32(n)
    var A = [Float](repeating: 0, count: n * n)
    memcpy(&A, matrix, n * n * 4)
    var LDA: Int32 = Int32(n)
    var W = [Float](repeating: 0, count: n)
    var WORK: [Float] = [0]
    var LWORK: Int32 = -1
    var IWORK: [Int32] = [0]
    var LIWORK: Int32 = -1
    var INFO: Int32 = 0
    A.withContiguousMutableStorageIfAvailable {
      let A = $0.baseAddress!
      W.withContiguousMutableStorageIfAvailable {
        let W = $0.baseAddress!
        WORK.withContiguousMutableStorageIfAvailable {
          let WORK = $0.baseAddress!
          IWORK.withContiguousMutableStorageIfAvailable {
            let IWORK = $0.baseAddress!
            ssyevd_(
              &JOBZ, // JOBZ
              &UPLO, // UPLO
              &N, // N
              A, // A
              &LDA, // LDA
              W, // W
              WORK, // WORK
              &LWORK, // LWORK
              IWORK, // IWORK
              &LIWORK, // LIWORK
              &INFO // INFO
            )
            guard INFO == 0 else {
              fatalError("LAPACK error code: \(INFO)")
            }
          }
        }
        
        LWORK = Int32(WORK[0])
        LIWORK = Int32(IWORK[0])
        WORK = [Float](repeating: 0, count: Int(LWORK))
        IWORK = [Int32](repeating: 0, count: Int(LIWORK))
        WORK.withContiguousMutableStorageIfAvailable {
          let WORK = $0.baseAddress!
          IWORK.withContiguousMutableStorageIfAvailable {
            let IWORK = $0.baseAddress!
            ssyevd_(
              &JOBZ, // JOBZ
              &UPLO, // UPLO
              &N, // N
              A, // A
              &LDA, // LDA
              W, // W
              WORK, // WORK
              &LWORK, // LWORK
              IWORK, // IWORK
              &LIWORK, // LIWORK
              &INFO // INFO
            )
            guard INFO == 0 else {
              fatalError("LAPACK error code: \(INFO)")
            }
          }
        }
      }
    }
    return (eigenvalues: W, eigenvectors: A)
  }
  
  // Solves a linear equation, where the matrix is in row-major order.
  static func solveLinearSystem(
    matrix: [Float], vector: [Float], n: Int
  ) -> [Float] {
    guard matrix.count == n * n,
          vector.count == n else {
      fatalError("Arguments to linear equation had invalid size.")
    }
    
    // Transpose the matrix into column-major order.
    var coefficients = [Float](repeating: .zero, count: n * n)
    for rowID in 0..<n {
      for columnID in 0..<n {
        let oldAddress = rowID * n + columnID
        let newAddress = columnID * n + rowID
        let value = matrix[oldAddress]
        coefficients[newAddress] = value
      }
    }
    
    // Call into the LAPACK function.
    var N: Int32 = Int32(n)
    var NRHS: Int32 = 1
    var A: [Float] = coefficients
    var LDA: Int32 = Int32(n)
    var IPIV: [Int32] = .init(repeating: .zero, count: n)
    var B: [Float] = vector
    var LDB: Int32 = Int32(n)
    var INFO: Int32 = 0
    A.withContiguousMutableStorageIfAvailable {
      let A = $0.baseAddress!
      B.withContiguousMutableStorageIfAvailable {
        let B = $0.baseAddress!
        sgesv_(
          &N,
          &NRHS,
          A,
          &LDA,
          &IPIV,
          B,
          &LDB,
          &INFO)
        guard INFO == 0 else {
          fatalError("LAPACK error code: \(INFO)")
        }
      }
    }
    
    // Return B, which was overwritten with the solution.
    let X = B
    return X
  }
}
