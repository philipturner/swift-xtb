import XCTest
import Accelerate
import Mechanosynthesis
import Numerics

final class DiagonalizationTests: XCTestCase {
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
    eigenvectors = LinearAlgebraUtilities
      .transpose(matrix: eigenvectors, n: 7)
    eigenvectors = LinearAlgebraUtilities
      .modifiedGramSchmidt(matrix: eigenvectors, n: 7)
    
    // Well-conditioned eigenspectra without degenerate clusters.
    testEigenvalues([
      4, 3, 2, 1, 0.1, -1.3, -2.3
    ], degenerateEigenvalues: [])
    testEigenvalues([
      4, 3, 2, 1, 0, -1, -2
    ], degenerateEigenvalues: [])
    
    // Problematic for the iterative diagonalization experiment.
    testEigenvalues([
      4, 3.01, 3.00, 2.99, 0, -2, -2
    ], degenerateEigenvalues: [1, 0])
    testEigenvalues([
      4, 3.5, 3, 2, 0, -1, -2
    ], degenerateEigenvalues: [])
    testEigenvalues([
      2, 2, 0, -2.99, -3.00, -3.01, -4
    ], degenerateEigenvalues: [6, 5])
    testEigenvalues([
      2, 1, 0, -2.99, -3.00, -3.01, -4
    ], degenerateEigenvalues: [])
    
    // Not problematic for the iterative diagonalization experiment.
    testEigenvalues([
      4, 3.01, 3.00, 2.99, 0, -1, -2
    ], degenerateEigenvalues: [])
    testEigenvalues([
      4, 3.5, 3, 2, 0, -2, -2
    ], degenerateEigenvalues: [1, 0])
    testEigenvalues([
      2, 2, 0, -2, -3, -3.5, -4
    ], degenerateEigenvalues: [6, 5])
    
    func testEigenvalues(
      _ eigenvalues: [Float],
      degenerateEigenvalues: Set<Int>
    ) {
      var Λ = [Float](repeating: 0, count: 7 * 7)
      for i in 0..<7 {
        let address = i * 7 + i
        let value = eigenvalues[i]
        Λ[address] = value
      }
      let ΣT = eigenvectors
      
      let ΛΣT = LinearAlgebraUtilities.matrixMultiply(
        matrixA: Λ, transposeA: false,
        matrixB: ΣT, transposeB: false, n: 7)
      let A = LinearAlgebraUtilities.matrixMultiply(
        matrixA: ΣT, transposeA: true,
        matrixB: ΛΣT, transposeB: false, n: 7)
      let AΣT = LinearAlgebraUtilities.matrixMultiply(
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
      
      // Check that Diagonalization properly handles this matrix.
      executeLAPACKComparison(
        matrix: A, n: 7, degenerateEigenvalues: degenerateEigenvalues)
    }
    
    func testMatrix(_ originalMatrixA: [Float], n: Int) {
      let originalMatrixT = LinearAlgebraUtilities
        .tridiagonalize(matrix: originalMatrixA, n: n)
      let (D, Z) = LinearAlgebraUtilities
        .divideAndConquer(matrix: originalMatrixT, n: n)
      
      // Check that the eigenvectors produce the eigenvalues.
      let HΨ = LinearAlgebraUtilities.matrixMultiply(
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
  
  // This test covers an algorithm for generating bulge chasing sequences.
  func testBulgeChasing() throws {
    testMatrix(n: 3, nb: 2)
    testMatrix(n: 4, nb: 2)
    testMatrix(n: 4, nb: 3)
    testMatrix(n: 6, nb: 2)
    testMatrix(n: 10, nb: 2)
    testMatrix(n: 10, nb: 3)
    testMatrix(n: 10, nb: 4)
    testMatrix(n: 10, nb: 8)
    testMatrix(n: 10, nb: 9)
    testMatrix(n: 11, nb: 2)
    testMatrix(n: 11, nb: 3)
    testMatrix(n: 11, nb: 4)
    testMatrix(n: 11, nb: 8)
    testMatrix(n: 11, nb: 9)
    testMatrix(n: 11, nb: 10)
    testMatrix(n: 19, nb: 5)
    testMatrix(n: 27, nb: 8)
    testMatrix(n: 32, nb: 8)
    testMatrix(n: 33, nb: 8)
    
    func testMatrix(n: Int, nb: Int) {
      var matrixA = [Int](repeating: 0, count: n * n)
      for vectorID in 0..<n {
        let address = vectorID * n + vectorID
        matrixA[address] = 1
        
        for subDiagonalID in 1...nb {
          if vectorID + subDiagonalID < n {
            let addressLower = (vectorID + subDiagonalID) * n + vectorID
            let addressUpper = vectorID * n + (vectorID + subDiagonalID)
            matrixA[addressUpper] = 1
            matrixA[addressLower] = 1
          }
        }
      }
      
      // [sweepID, startColumnID, startRowID, endRowID]
      var bulgeChasingSequence: [SIMD4<Int>] = []
      
      func startSweep(sweepID: Int) {
        let startVectorID = sweepID + 1
        guard startVectorID < n - 1 else {
          fatalError("Attempted a sweep that will not generate any bulges.")
        }
        let endVectorID = min(sweepID + nb + 1, n)
        bulgeChasingSequence.append(SIMD4(
          sweepID, sweepID, startVectorID, endVectorID))
        
        var nextBulgeCornerID: SIMD2<Int>?
        nextBulgeCornerID = applyBulgeChase(
          startColumnID: sweepID,
          startRowID: startVectorID,
          endRowID: endVectorID)
        while nextBulgeCornerID != nil {
          let cornerRowID = nextBulgeCornerID![0]
          let cornerColumnID = nextBulgeCornerID![1]
          let endVectorID = cornerRowID + 1
          let startVectorID = cornerColumnID + nb
          
          if endVectorID - startVectorID > 1 {
            bulgeChasingSequence.append(SIMD4(
              sweepID, cornerColumnID, startVectorID, endVectorID))
            
            nextBulgeCornerID = applyBulgeChase(
              startColumnID: cornerColumnID,
              startRowID: startVectorID,
              endRowID: endVectorID)
          } else {
            nextBulgeCornerID = nil
          }
        }
      }
      
      // Returns the location of the next bulge to chase.
      // - startColumnID: The reflector projected onto the main diagonal.
      // - startRowID: The element ID to pivot on.
      // - endRowID: One past the last element to affect.
      func applyBulgeChase(
        startColumnID: Int,
        startRowID: Int,
        endRowID: Int
      ) -> SIMD2<Int>? {
        // Apply Householder reflections to the first column.
        for rowID in startRowID..<endRowID {
          if rowID != startRowID {
            let addressLower = rowID * n + startColumnID
            let addressUpper = startColumnID * n + rowID
            matrixA[addressLower] = 0
            matrixA[addressUpper] = 0
          }
        }
        
        // Loop until you don't find any more ones.
        var nextBulgeCornerID: SIMD2<Int>? // (row, column)
        for columnID in (startColumnID + 1)..<n {
          var foundOne = false
          for rowID in startRowID..<endRowID {
            let address = rowID * n + columnID
            if matrixA[address] == 1 {
              foundOne = true
            }
          }
          guard foundOne else {
            break
          }
          for rowID in startRowID..<endRowID {
            let addressLower = rowID * n + columnID
            let addressUpper = columnID * n + rowID
            matrixA[addressLower] = 1
            matrixA[addressUpper] = 1
          }
          
          // The column and row ID refer to the coordinates on the other side
          // of the symmetric matrix.
          nextBulgeCornerID = SIMD2(columnID, startRowID)
        }
        return nextBulgeCornerID
      }
      
      guard nb > 1 else {
        // The above chaser failed for this edge case. Bake the zero-operation
        // condition into the algorithm for predicting sequences.
        fatalError("Cannot perform bulge chasing when nb = 1.")
      }
      for sweepID in 0..<max(0, n - 2) {
        startSweep(sweepID: sweepID)
      }
      
      var chasingOperationCursor = 0
      func assertNextOperation(
        sweepID: Int, startColumnID: Int,
        startRowID: Int, endRowID: Int
      ) {
        guard chasingOperationCursor < bulgeChasingSequence.count else {
          XCTFail("Overflowed the chasing sequence buffer.")
          return
        }
        let currentOperation = bulgeChasingSequence[chasingOperationCursor]
        
        XCTAssertEqual(
          sweepID, currentOperation[0],
          "sequence[\(chasingOperationCursor)]/sweepID")
        XCTAssertEqual(
          startColumnID, currentOperation[1],
          "sequence[\(chasingOperationCursor)]/startColumnID")
        XCTAssertEqual(
          startRowID, currentOperation[2],
          "sequence[\(chasingOperationCursor)]/startRowID")
        XCTAssertEqual(
          endRowID, currentOperation[3],
          "sequence[\(chasingOperationCursor)]/endRowID")
        
        chasingOperationCursor += 1
      }
      
      for sweepID in 0..<max(0, n - 2) {
        let startVectorID = sweepID + 1
        var endVectorID = sweepID + nb + 1
        endVectorID = min(endVectorID, n)
        guard endVectorID - startVectorID > 1 else {
          fatalError("Generated empty Householder transform.")
        }
        assertNextOperation(
          sweepID: sweepID, startColumnID: sweepID,
          startRowID: startVectorID, endRowID: endVectorID)
        
        var operationID = 1
        while true {
          let startColumnID = (sweepID - nb + 1) + operationID * nb
          let startVectorID = (sweepID + 1) + operationID * nb
          var endVectorID = (sweepID + nb + 1) + operationID * nb
          endVectorID = min(endVectorID, n)
          
          if endVectorID - startVectorID > 1 {
            assertNextOperation(
              sweepID: sweepID, startColumnID: startColumnID,
              startRowID: startVectorID, endRowID: endVectorID)
          } else {
            break
          }
          operationID += 1
        }
      }
      XCTAssertEqual(chasingOperationCursor, bulgeChasingSequence.count)
    }
  }
  
  // Test the two-stage process for tridiagonalizing a matrix, and the
  // aggregation of the Householder reflectors.
  func testTwoStageTridiagonalization() throws {
    testBlockSize(nb: 1, sb: 1)
    testBlockSize(nb: 2, sb: 1)
    testBlockSize(nb: 3, sb: 1)
    testBlockSize(nb: 4, sb: 1)
    testBlockSize(nb: 5, sb: 1)
    testBlockSize(nb: 6, sb: 1)
    
    testBlockSize(nb: 2, sb: 2)
    testBlockSize(nb: 3, sb: 3)
    testBlockSize(nb: 4, sb: 4)
    testBlockSize(nb: 5, sb: 5)
    testBlockSize(nb: 6, sb: 6)
    
    testBlockSize(nb: 4, sb: 2)
    testBlockSize(nb: 6, sb: 2)
    testBlockSize(nb: 6, sb: 3)
    
    func testBlockSize(nb: Int, sb: Int) {
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
      
      // Make the matrix symmetric.
      originalMatrixA = LinearAlgebraUtilities.matrixMultiply(
        matrixA: originalMatrixA,
        matrixB: originalMatrixA,
        transposeB: true, n: n)
      
      var diagonalizationDesc = DiagonalizationDescriptor()
      diagonalizationDesc.matrix = originalMatrixA
      diagonalizationDesc.problemSize = n
      diagonalizationDesc.blockSize = nb
      diagonalizationDesc.smallBlockSize = sb
      
      let diagonalization = Diagonalization(descriptor: diagonalizationDesc)
      let eigenvalues = diagonalization.eigenvalues
      let eigenvectors = diagonalization.eigenvectors
      
      let expectedEigenvalues: [Float] = [
        0.0011429853, 0.5075689, 0.75631386, 6.188949, 145.55783, 443.30386,
        7871.384
      ]
      let eigenvalueAccuracies: [Float] = [
        1e-4, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-2
      ]
      
      var expectedEigenvectors: [[Float]] = []
      expectedEigenvectors.append([
        0.0012315774, 0.09477102, 0.08277881, -0.8589063, -0.49401814,
        0.047964625, -0.0094997585])
      expectedEigenvectors.append([
        -0.40089795, 0.4541109, 0.13252003, 0.43725404, -0.640054,
         0.120887995, -0.0053356886])
      expectedEigenvectors.append([
        -0.55025303, 0.54824287, 0.056886543, -0.26165208, 0.56977874,
         0.016768508, 0.004901767])
      expectedEigenvectors.append([
        0.30370802, 0.07415314, 0.7463128, 0.029480757, 0.14456089,
        0.5677579, -0.034115434])
      expectedEigenvectors.append([
        0.5517132, 0.57923084, -0.51686704, 0.016868684, 0.026573017,
        0.2974512, -0.059032507])
      expectedEigenvectors.append([
        -0.3680466, -0.37356234, -0.38447195, -0.038388584, 0.0019791797,
         0.75624263, 0.06159757])
      expectedEigenvectors.append([
        0.06645671, 0.06063038, 0.019931633, -0.00017843395, -0.0045418143,
        -0.008672804, 0.99569774])
      
      for i in 0..<7 {
        let expected = expectedEigenvalues[i]
        let actual = eigenvalues[i]
        let accuracy = eigenvalueAccuracies[i]
        XCTAssertEqual(actual, expected, accuracy: accuracy)
      }
      
      for i in 0..<7 {
        let expected = expectedEigenvectors[i]
        var actual = [Float](repeating: 0, count: 7)
        for elementID in 0..<7 {
          let address = i * 7 + elementID
          actual[elementID] = eigenvectors[address]
        }
        
        var dotProduct: Float = .zero
        for elementID in 0..<7 {
          dotProduct += expected[elementID] * actual[elementID]
        }
        XCTAssertEqual(dotProduct.magnitude, 1, accuracy: 1e-5)
      }
    }
  }
  
  // Test every edge case that relates to small matrices.
  func testSmallMatrixDiagonalization() throws {
    // n = 1
    executeLAPACKComparison(matrix: [1], n: 1)
    executeLAPACKComparison(matrix: [2], n: 1)
    executeLAPACKComparison(matrix: [3], n: 1)
    executeLAPACKComparison(matrix: [0], n: 1)
    executeLAPACKComparison(matrix: [0.0001], n: 1)
    executeLAPACKComparison(matrix: [-1], n: 1)
    
    // n = 2
    executeLAPACKComparison(matrix: [
      1, -4,
      -4, 9,
    ], n: 2)
    executeLAPACKComparison(matrix: [
      1, 4,
      4, 9,
    ], n: 2)
    executeLAPACKComparison(matrix: [
      -1, -4,
      -4, -9,
    ], n: 2)
    executeLAPACKComparison(matrix: [
      0, 2,
      2, 3,
    ], n: 2)
    executeLAPACKComparison(matrix: [
      0, 0,
      0, 1,
    ], n: 2)
    executeLAPACKComparison(matrix: [
      0, 1,
      1, 0,
    ], n: 2)
    executeLAPACKComparison(matrix: [
      0, 0,
      0, 0,
    ], n: 2, degenerateEigenvalues: [0, 1])
    executeLAPACKComparison(matrix: [
      0, 0,
      0, 0.0001,
    ], n: 2)
    executeLAPACKComparison(matrix: [
      0, -10,
      -10, 0.0001,
    ], n: 2)
    
    // n = 3
    executeLAPACKComparison(matrix: [
      0, 0, 1,
      0, 2, 0,
      1, 0, 0,
    ], n: 3)
    executeLAPACKComparison(matrix: [
      8, 7, 6,
      7, -2, -3,
      6, -3, -3,
    ], n: 3)
    executeLAPACKComparison(matrix: [
      8, 0, 6,
      0, -2, 0,
      6, 0, -3,
    ], n: 3)
    executeLAPACKComparison(matrix: [
      8, -0, 6,
      -0, -2, -0,
      6, -0, -3,
    ], n: 3)
    
    // n = 4
    executeLAPACKComparison(matrix: [
      4, 0.1, 0, 0,
      0.1, 3, 0.1, 0,
      0, 0.1, 2, 0.1,
      0, 0, 0.1, 1,
    ], n: 4)
    executeLAPACKComparison(matrix: [
      4, 0, 0, 0,
      0, 3, 0, 0,
      0, 0, 2, 0,
      0, 0, 0, 1,
    ], n: 4)
    executeLAPACKComparison(matrix: [
      1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 1, 0,
      0, 0, 0, 1,
    ], n: 4, degenerateEigenvalues: [0, 1, 2, 3])
  }
  
  // Test 2nd-order and 4th-order finite difference operators.
  func testLaplacian() throws {
    // Laplacian matrix (4x4)
    executeLAPACKComparison(matrix: [
      2, -1, 0, 0,
      -1, 2, -1, 0,
      0, -1, 2, -1,
      0, 0, -1, 2,
    ], n: 4)
    executeLAPACKComparison(matrix: [
      -2, 1, 0, 1,
      1, -2, 1, 0,
      0, 1, -2, 1,
      1, 0, 1, -2,
    ], n: 4, degenerateEigenvalues: [1, 2])
    executeLAPACKComparison(matrix: [
      2, -1, 0, -1,
      -1, 2, -1, 0,
      0, -1, 2, -1,
      -1, 0, -1, 2,
    ], n: 4, degenerateEigenvalues: [1, 2])
    
    // Laplacian matrix (5x5)
    executeLAPACKComparison(matrix: [
      -2, 1, 0, 0, 0,
      1, -2, 1, 0, 0,
      0, 1, -2, 1, 0,
      0, 0, 1, -2, 1,
      0, 0, 0, 1, -2,
    ], n: 5)
    executeLAPACKComparison(matrix: [
      -2, 1, 0, 0, 1,
      1, -2, 1, 0, 0,
      0, 1, -2, 1, 0,
      0, 0, 1, -2, 1,
      1, 0, 0, 1, -2,
    ], n: 5, degenerateEigenvalues: [0, 1, 2, 3])
    
    // Laplacian matrix (6x6)
    executeLAPACKComparison(matrix: [
      -2, 1, 0, 0, 0, 0,
      1, -2, 1, 0, 0, 0,
      0, 1, -2, 1, 0, 0,
      0, 0, 1, -2, 1, 0,
      0, 0, 0, 1, -2, 1,
      0, 0, 0, 0, 1, -2,
    ], n: 6)
    executeLAPACKComparison(matrix: [
      -2, 1, 0, 0, 0, 1,
      1, -2, 1, 0, 0, 0,
      0, 1, -2, 1, 0, 0,
      0, 0, 1, -2, 1, 0,
      0, 0, 0, 1, -2, 1,
      1, 0, 0, 0, 1, -2,
    ], n: 6, degenerateEigenvalues: [1, 2, 3, 4])
    
    // Fourth-order finite differencing
    executeLAPACKComparison(matrix: [
      6, -4, 1, 0, 0, 0,
      -4, 6, -4, 1, 0, 0,
      1, -4, 6, -4, 1, 0,
      0, 1, -4, 6, -4, 1,
      0, 0, 1, -4, 6, -4,
      0, 0, 0, 1, -4, 6,
    ], n: 6)
    executeLAPACKComparison(matrix: [
      6, -4, 1, 0, 1, -4,
      -4, 6, -4, 1, 0, 1,
      1, -4, 6, -4, 1, 0,
      0, 1, -4, 6, -4, 1,
      1, 0, 1, -4, 6, -4,
      -4, 1, 0, 1, -4, 6,
    ], n: 6, degenerateEigenvalues: [1, 2, 3, 4])
  }
  
  // Reproduce the LAPACK benchmarking experiment where Accelerate's
  // two-stage solver failed.
  func testBenchmark() throws {
    func benchmarkProblemSize(n: Int, trialCount: Int) {
      // Generate the eigenvalues.
      var Λ = [Float](repeating: 0, count: n * n)
      for vectorID in 0..<n {
        let position = Float(n / 2 - vectorID)
        var sign: Float
        if position > 0 {
          sign = 1
        } else if position == 0 {
          sign = 0
        } else {
          sign = -1
        }
        
        let powerPart = Float.pow(1.01, position.magnitude)
        let address = vectorID * n + vectorID
        Λ[address] = sign * powerPart
      }
      
      // Generate the eigenvectors.
      var Σ = [Float](repeating: 0, count: n * n)
      for vectorID in 0..<n {
        var vector = [Float](repeating: 0, count: n)
        for elementID in 0..<n {
          vector[elementID] = .random(in: -1...1)
        }
        for elementID in 0..<n {
          let address = elementID * n + vectorID
          Σ[address] = vector[elementID]
        }
      }
      Σ = LinearAlgebraUtilities
        .modifiedGramSchmidt(matrix: Σ, n: n)
      
      // Construct the symmetric matrix.
      let ΛΣT = LinearAlgebraUtilities.matrixMultiply(
        matrixA: Λ, matrixB: Σ, transposeB: true, n: n)
      let A = LinearAlgebraUtilities.matrixMultiply(
        matrixA: Σ, matrixB: ΛΣT, n: n)
      let AΣ = LinearAlgebraUtilities.matrixMultiply(
        matrixA: A, matrixB: Σ, n: n)
      
      // Check self-consistency of the matrix, prior to diagonalization.
      var eigenvaluePassRate: Float = .zero
      for vectorID in 0..<n {
        let expected = Λ[vectorID * n + vectorID]
        var actual: Float = .zero
        for elementID in 0..<n {
          // Accumulate the Rayleigh quotient to find the eigenvalue.
          let address = elementID * n + vectorID
          actual += Σ[address] * AΣ[address]
        }
        
        let accuracy = max(1e-5, expected.magnitude * 1e-3)
        let error = (actual - expected).magnitude
        if error.magnitude < accuracy {
          eigenvaluePassRate += 1 / Float(n)
        }
      }
      XCTAssertGreaterThan(eigenvaluePassRate, 0.90)
      
      var diagonalizationDesc = DiagonalizationDescriptor()
      diagonalizationDesc.matrix = A
      diagonalizationDesc.problemSize = n
      if n < 8 {
        diagonalizationDesc.blockSize = n
        diagonalizationDesc.smallBlockSize = 1
      } else if n < 32 {
        diagonalizationDesc.blockSize = 8
        diagonalizationDesc.smallBlockSize = 2
      } else {
        diagonalizationDesc.blockSize = 32
        diagonalizationDesc.smallBlockSize = 8
      }
      let diagonalization = Diagonalization(descriptor: diagonalizationDesc)
      
      let oneStageEigenvalues = LinearAlgebraUtilities
        .diagonalize(matrix: A, n: n).0
      #if false
      let twoStageEigenvalues = LinearAlgebraUtilities
        .eigenvaluesTwoStage(matrix: A, n: n)
      #else
      let twoStageEigenvalues = diagonalization.eigenvalues
      #endif
      
      for i in 0..<n {
        let expectedID = (n - 1) - i
        let expected = Λ[expectedID * n + expectedID]
        let actual1 = oneStageEigenvalues[i]
        let actual2 = twoStageEigenvalues[i]
        let accuracy = max(1e-5, expected.magnitude * 1e-3)
        XCTAssertEqual(
          actual1, expected, accuracy: accuracy,
          "One-stage solver failed at Λ[\(i)].")
        XCTAssertEqual(
          actual2, expected, accuracy: accuracy,
          "Two-stage solver failed at Λ[\(i)].")
        XCTAssertEqual(
          actual1, actual2, accuracy: accuracy,
          "Solvers disagreed at Λ[\(i)].")
      }
    }
    
    // A small test that the benchmark works correctly.
    benchmarkProblemSize(n: 7, trialCount: 1)
    for n in 1...20 {
      benchmarkProblemSize(n: n, trialCount: 3)
    }
    for nTenth in 3...10 {
      benchmarkProblemSize(n: 10 * nTenth, trialCount: 3)
    }
  }
}

// Diagonalizes the matrix and checks that the results agree with LAPACK.
// - degenerateEigenvalues: the indices of degenerate eigenpairs, in
//                          ascending order of the signed eigenvalues.
private func executeLAPACKComparison(
  matrix: [Float],
  n: Int,
  degenerateEigenvalues: Set<Int> = []
) {
  let (expectedEigenvalues, expectedEigenvectors) = LinearAlgebraUtilities
    .diagonalize(matrix: matrix, n: n)
  
  var diagonalizationDesc = DiagonalizationDescriptor()
  diagonalizationDesc.matrix = matrix
  diagonalizationDesc.problemSize = n
  if n <= 4 {
    diagonalizationDesc.blockSize = 2
  } else {
    diagonalizationDesc.blockSize = 4
  }
  
  let diagonalizaton = Diagonalization(descriptor: diagonalizationDesc)
  let actualEigenvalues = diagonalizaton.eigenvalues
  let actualEigenvectors = diagonalizaton.eigenvectors
  
  for i in 0..<n {
    let expected = expectedEigenvalues[i]
    let actual = actualEigenvalues[i]
    let accuracy = max(1e-3 * expected.magnitude, 1e-5)
    XCTAssertEqual(actual, expected, accuracy: accuracy)
  }
  for i in 0..<n {
    if degenerateEigenvalues.contains(i) {
      continue
    }
    
    var expected = [Float](repeating: 0, count: n)
    var actual = [Float](repeating: 0, count: n)
    for elementID in 0..<n {
      let address = i * n + elementID
      expected[elementID] = expectedEigenvectors[address]
      actual[elementID] = actualEigenvectors[address]
    }
    
    var dotProduct: Float = .zero
    for elementID in 0..<n {
      dotProduct += expected[elementID] * actual[elementID]
    }
    XCTAssertEqual(dotProduct.magnitude, 1, accuracy: 1e-5)
  }
}
