//
//  DivideAndConquer.swift
//
//
//  Created by Philip Turner on 3/11/24.
//

import Accelerate

extension Diagonalization {
  // Diagonalizes a tridiagonal matrix with LAPACK divide-and-conquer.
  // Returns the eigenvectors in column-major format, as output by LAPACK.
  mutating func solveTridiagonalEigenproblem() {
    // Store the tridiagonal matrix in a compact form.
    var D = [Float](repeating: 0, count: problemSize)
    var E = [Float](repeating: 0, count: problemSize - 1)
    for diagonalID in 0..<problemSize {
      let matrixAddress = diagonalID * problemSize + diagonalID
      let vectorAddress = diagonalID
      D[vectorAddress] = matrix[matrixAddress]
    }
    for subDiagonalID in 0..<problemSize - 1 {
      let rowID = subDiagonalID
      let columnID = subDiagonalID + 1
      let matrixAddress = rowID * problemSize + columnID
      let vectorAddress = rowID
      E[vectorAddress] = matrix[matrixAddress]
    }
    
    // Query the workspace size.
    var JOBZ = CChar(Character("V").asciiValue!)
    var N = Int32(problemSize)
    var LDZ = Int32(problemSize)
    var WORK = [Float](repeating: 0, count: 1)
    var LWORK = Int32(-1)
    var IWORK = [Int32](repeating: 0, count: 1)
    var LIWORK = Int32(-1)
    var INFO = Int32(0)
    sstevd_(
      &JOBZ, &N, nil, nil, nil, &LDZ, &WORK, &LWORK, &IWORK, &LIWORK, &INFO)
    
    // Call into LAPACK.
    var Z = [Float](repeating: 0, count: problemSize * problemSize)
    LWORK = Int32(WORK[0])
    LIWORK = Int32(IWORK[0])
    WORK = [Float](repeating: 0, count: Int(LWORK))
    IWORK = [Int32](repeating: 0, count: Int(LIWORK))
    sstevd_(
      &JOBZ, &N, &D, &E, &Z, &LDZ, &WORK, &LWORK, &IWORK, &LIWORK, &INFO)
    
    // Return the eigenpairs.
    eigenvalues = D
    eigenvectors = Z
  }
}
