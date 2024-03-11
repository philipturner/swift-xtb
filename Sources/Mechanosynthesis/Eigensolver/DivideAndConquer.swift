//
//  DivideAndConquer.swift
//
//
//  Created by Philip Turner on 3/11/24.
//

import Accelerate

// Diagonalizes a tridiagonal matrix with LAPACK divide-and-conquer.
// Returns the eigenvectors in column-major format, as output by LAPACK.
func divideAndConquer(matrix: [Float], n: Int) -> (
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
