//
//  main.swift
//
//
//  Created by Philip Turner on 5/25/24.
//

import Accelerate
import Foundation
import LinearAlgebra
import Numerics
import xTB

// Set up a test for various matrices with different contents and sizes.
// - How to fill the matrices (analytical eigenpairs)
//   - Matrix multiply, transpose, and modified GS utility functions
// - How to check accuracy of matrices
//   - Copy algorithm from test suite for suitability of the test, prior to
//     actual testing
//   - Check the eigenvalues only (not the eigenvectors)
// - After getting the correctness test, we can transform it into a performance
//   test. Finally, a performance test that tests multiple problem sizes and
//   formats the data.

// Setup procedure: eigenvalues on a logarithmic spectrum. We'll need to modify
// this to also process the identity matrix.
enum TestType {
  case identityMatrix
  case logarithmicSpectrum
}

struct TestDescriptor {
  var problemSize: Int?
  var type: TestType?
}

struct TestCase {
  // A square matrix of eigenvalues, NOT a linear array of eigenvalues.
  var eigenvalues: [Float]
  
  // A square matrix of eigenvectors.
  var eigenvectors: [Float]
  
  // A hamiltonian matrix.
  var matrix: [Float]
  
  // The side length of the square, in array slots.
  var problemSize: Int
  
  init(descriptor: TestDescriptor) {
    guard let problemSize = descriptor.problemSize,
          let type = descriptor.type else {
      fatalError("Descriptor was incomplete.")
    }
    self.problemSize = problemSize
    
    switch type {
    case .identityMatrix:
      eigenvalues = [Float](repeating: 0, count: problemSize * problemSize)
      eigenvectors = [Float](repeating: 0, count: problemSize * problemSize)
      
      for diagonalID in 0..<problemSize {
        let address = diagonalID * problemSize + diagonalID
        eigenvalues[address] = 1
        eigenvectors[address] = 1
      }
      
    case .logarithmicSpectrum:
      // Generate the eigenvalues.
      eigenvalues = [Float](repeating: 0, count: problemSize * problemSize)
      for vectorID in 0..<problemSize {
        let position = Float(problemSize / 2 - vectorID)
        var sign: Float
        if position > 0 {
          sign = 1
        } else if position == 0 {
          sign = 0
        } else {
          sign = -1
        }
        
        let powerPart = Float.pow(1.01, position.magnitude)
        let address = vectorID * problemSize + vectorID
        eigenvalues[address] = sign * powerPart
      }
      
      // Generate the eigenvectors.
      eigenvectors = [Float](repeating: 0, count: problemSize * problemSize)
      for vectorID in 0..<problemSize {
        var vector = [Float](repeating: 0, count: problemSize)
        for elementID in 0..<problemSize {
          vector[elementID] = .random(in: -1...1)
        }
        for elementID in 0..<problemSize {
          let address = elementID * problemSize + vectorID
          eigenvectors[address] = vector[elementID]
        }
      }
      eigenvectors = LinearAlgebraUtilities
        .modifiedGramSchmidt(matrix: eigenvectors, n: problemSize)
    }
    
    // Construct the symmetric matrix.
    let ΛΣT = LinearAlgebraUtilities.matrixMultiply(
      matrixA: eigenvalues, matrixB: eigenvectors,
      transposeB: true, n: problemSize)
    matrix = LinearAlgebraUtilities.matrixMultiply(
      matrixA: eigenvectors, matrixB: ΛΣT, n: problemSize)
    
    // Check that you just generated a valid matrix.
    checkSelfConsistency()
  }
  
  // Verify that a test will actually accept the eigenproblem.
  private func checkSelfConsistency() {
    let AΣ = LinearAlgebraUtilities.matrixMultiply(
      matrixA: matrix, matrixB: eigenvectors, n: problemSize)
    
    var eigenvaluePassRate: Float = .zero
    for vectorID in 0..<problemSize {
      let expected = eigenvalues[vectorID * problemSize + vectorID]
      var actual: Float = .zero
      for elementID in 0..<problemSize {
        // Accumulate the Rayleigh quotient to find the eigenvalue.
        let address = elementID * problemSize + vectorID
        actual += eigenvectors[address] * AΣ[address]
      }
      
      let accuracy = max(1e-5, expected.magnitude * 1e-3)
      let error = (actual - expected).magnitude
      if error.magnitude < accuracy {
        eigenvaluePassRate += 1 / Float(problemSize)
      }
    }
    guard eigenvaluePassRate > 0.90 else {
      fatalError("Problem is not working correctly.")
    }
  }
  
  // Verify the eigenvalues against a different answer.
  //
  // NOTE: The entered eigenvalues must be a linear array. This differs from
  // the memory format that the test uses (a square matrix).
  func check(eigenvalues otherEigenvalues: [Float]) {
    for i in 0..<problemSize {
      let expectedID = (problemSize - 1) - i
      let expected = eigenvalues[expectedID * problemSize + expectedID]
      let actual = otherEigenvalues[i]
      let accuracy = max(1e-5, expected.magnitude * 1e-3)
      guard (actual - expected).magnitude < accuracy else {
        fatalError("Solver failed at  Λ[\(i)].")
      }
    }
  }
}

// Double precision implementation of diagonalization.
func diagonalize(matrix: [Double], n: Int, jobz: Character) -> (
  eigenvalues: [Double], eigenvectors: [Double]
) {
  var JOBZ = CChar(jobz.asciiValue!)
  var UPLO = CChar(Character("L").asciiValue!)
  var N: Int32 = Int32(n)
  var A = [Double](repeating: 0, count: n * n)
  memcpy(&A, matrix, n * n * 8)
  var LDA: Int32 = Int32(n)
  var W = [Double](repeating: 0, count: n)
  var WORK: [Double] = [0]
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
          dsyevd_(
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
      WORK = [Double](repeating: 0, count: Int(LWORK))
      IWORK = [Int32](repeating: 0, count: Int(LIWORK))
      WORK.withContiguousMutableStorageIfAvailable {
        let WORK = $0.baseAddress!
        IWORK.withContiguousMutableStorageIfAvailable {
          let IWORK = $0.baseAddress!
          dsyevd_(
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

// Generate the problem sizes to target.
var problemSizesRaw: [Int] = [
  147, 155, 157, 165, 167, 175, 177, 185, 187, 193, 201, 203, 211, 213, 219,
  221, 221, 222, 228, 230, 228, 235, 237, 244, 246, 253, 261, 263, 265, 265,
  271, 273, 279, 281, 287, 289, 297, 299, 306, 308, 316, 317, 325, 327
]
var problemSizes: [Int]
do {
  var problemSizesSet = Set(problemSizesRaw)
  for problemSizeID in problemSizesRaw.indices {
    var problemSize = problemSizesRaw[problemSizeID]
    problemSize += 80
    problemSizesSet.insert(problemSize)
  }
  problemSizes = Array(problemSizesSet).sorted()
}

// Run the benchmarks.
for problemSize in problemSizes {
  var testDescriptor = TestDescriptor()
  testDescriptor.problemSize = problemSize
  testDescriptor.type = .logarithmicSpectrum
  
  let test = TestCase(descriptor: testDescriptor)
  let matrix64 = test.matrix.map(Double.init)
  var eigenvalues64: [Double] = []
  
  var maxGFLOPS: Float = .zero
  for _ in 0..<3 {
    let checkpoint0 = Date()
    (eigenvalues64, _) = diagonalize(
      matrix: matrix64, n: problemSize, jobz: "N")
    let checkpoint1 = Date()
    
    let latency = checkpoint1.timeIntervalSince(checkpoint0) as Double
    let operationCount = problemSize * problemSize * problemSize
    let gflops = Float(operationCount) / Float(latency) / 1e9
    maxGFLOPS = max(maxGFLOPS, gflops)
  }
  
  // TODO: Round to one significant digit.
  print("problemSize = \(problemSize), maxGFLOPS = \(maxGFLOPS)")
  
  let eigenvalues = eigenvalues64.map(Float.init)
  test.check(eigenvalues: eigenvalues)
}
