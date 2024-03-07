import XCTest
import Accelerate // Gate out once there are Swift kernels for CPUs without AMX.
import Numerics
import QuartzCore

final class LinearAlgebraTests: XCTestCase {
  // MARK: - Linear Algebra Functions
  
  // Multiplies two square matrices.
  // - Mixed precision instead of single precision.
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
  
  // Forms an orthogonal basis of the matrix's rows.
  // - Mixed precision instead of single precision.
  // - Classical GS instead of modified GS.
  static func rowGramSchmidt(
    matrix originalMatrix: [Float], n: Int
  ) -> [Float] {
    // Operate on the (output) matrix in-place.
    var matrix = originalMatrix
    
    func normalize(electronID: Int) {
      var norm: Double = .zero
      for cellID in 0..<n {
        let value = matrix[electronID * n + cellID]
        norm += Double(value * value)
      }
      
      let normalizationFactor = 1 / norm.squareRoot()
      for cellID in 0..<n {
        var value = Double(matrix[electronID * n + cellID])
        value *= normalizationFactor
        matrix[electronID * n + cellID] = Float(value)
      }
    }
    
    for electronID in 0..<n {
      // Normalize the vectors before taking dot products.
      normalize(electronID: electronID)
    }
    
    for electronID in 0..<n {
      // Determine the magnitude of components parallel to previous vectors.
      var dotProducts = [Double](repeating: .zero, count: n)
      for neighborID in 0..<electronID {
        var dotProduct: Double = .zero
        for cellID in 0..<n {
          let value1 = matrix[electronID * n + cellID]
          let value2 = matrix[neighborID * n + cellID]
          dotProduct += Double(value1) * Double(value2)
        }
        dotProducts[neighborID] = dotProduct
      }
      
      // Subtract all components parallel to previous vectors.
      for cellID in 0..<n {
        var value1 = Double(matrix[electronID * n + cellID])
        for neighborID in 0..<electronID {
          let dotProduct = dotProducts[neighborID]
          let value2 = matrix[neighborID * n + cellID]
          value1 -= dotProduct * Double(value2)
        }
        matrix[electronID * n + cellID] = Float(value1)
      }
      
      // Rescale the orthogonal component to unit vector length.
      normalize(electronID: electronID)
    }
    
    return matrix
  }
  
  // Standalone function for tridiagonalizing a matrix.
  // - Performs several loop iterations with O(n^2) complexity each.
  // - Overall, the procedure has O(n^3) complexity.
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
  
  // Test the QR algorithm for finding eigenvalues of a tridiagonal matrix.
  func testQRAlgorithm() throws {
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
    
    // Predict what the Q^T A Q / Q^T T Q algorithm should produce.
    var expectedQ1 = originalMatrixT
    var expectedQ1T = Self.transpose(matrix: expectedQ1, n: n)
    expectedQ1T = Self.rowGramSchmidt(matrix: expectedQ1T, n: n)
    expectedQ1 = Self.transpose(matrix: expectedQ1T, n: n)
    print()
    print("expected Q1")
    for rowID in 0..<n {
      for columnID in 0..<n {
        let address = rowID * n + columnID
        let value = expectedQ1[address]
        print(value, terminator: ", ")
      }
      print()
    }
    
    var expectedT2 = Self.matrixMultiply(
      matrixA: originalMatrixT, matrixB: expectedQ1, n: n)
    expectedT2 = Self.matrixMultiply(
      matrixA: expectedQ1T, matrixB: expectedT2, n: n)
    print()
    print("expected T2")
    for rowID in 0..<n {
      for columnID in 0..<n {
        let address = rowID * n + columnID
        let value = expectedT2[address]
        print(value, terminator: ", ")
      }
      print()
    }
    
    // Store the tridiagonal matrix in a compact form.
    let currentMatrixT = originalMatrixT
    var a = [Float](repeating: 0, count: n + 1)
    var bSquared = [Float](repeating: 0, count: n)
    for diagonalID in 0..<n {
      let matrixAddress = diagonalID * n + diagonalID
      let vectorAddress = diagonalID
      a[vectorAddress] = currentMatrixT[matrixAddress]
    }
    for subDiagonalID in 0..<n - 1 {
      let rowID = subDiagonalID
      let columnID = subDiagonalID + 1
      let matrixAddress = rowID * n + columnID
      let vectorAddress = rowID
      
      let bValue = currentMatrixT[matrixAddress]
      bSquared[vectorAddress] = bValue * bValue
    }
    
    // Perform one iteration of QR decomposition.
    // https://doi.org/10.1093/comjnl/6.1.99
    for _ in 0..<5 {
      // Shift the matrix: T -> T - λI
      let λ = a[n - 1]
      for diagonalID in 0..<n {
        a[diagonalID] -= λ
      }
      
      // Factorize Q, R and compute RQ.
      print()
      print("iterating along diagonal")
      var u: Float = .zero
      var sSquared: Float = .zero
      var oneMinusSSquared: Float = 1
      for i in 0..<n {
        let γ = a[i] - u
        
        // Branch on the value of s^2.
//        var pSquared: Float
//        if (sSquared - 1).magnitude < 1e-6 {
//          fatalError("Edge case: s^2 converged to 1.")
//        } else {
//
//        }
        let pSquared = (γ * γ) / oneMinusSSquared
        
        // Branch on the value of i.
        let currentBSquared = bSquared[i]
        if i > 0 {
          bSquared[i - 1] = sSquared * (pSquared + currentBSquared)
        }
        sSquared = currentBSquared / (pSquared + currentBSquared)
        oneMinusSSquared = pSquared / (pSquared + currentBSquared)
        
        // Prepare for the next loop iteration.
        u = sSquared * (γ + a[i + 1])
        a[i] = γ + u
        
        // Log all of the state variables at the end of this iteration.
        print(i, "| γ = \(γ)")
        print(i, "| p^2 = \(pSquared)")
        print(i, "| s^2 = \(sSquared)")
        print(i, "| u = \(u)")
      }
      
      // Un-shift the matrix: RQ -> RQ + λI
      for diagonalID in 0..<n {
        a[diagonalID] += λ
      }
      
      print()
      print("actual a")
      for slotID in a.indices {
        print(slotID, "| a =", a[slotID])
      }
      
      print()
      print("actual b")
      for slotID in bSquared.indices {
        let value = bSquared[slotID]
        print(slotID, "| b^2 =", value, terminator: ", ")
        print("|b| =", value.squareRoot())
      }
    }
    
    /*
     Without shifting:
     
     iterating along diagonal
     0 | γ = 4.0
     0 | p^2 = 16.0
     0 | s^2 = 0.35999995
     0 | u = 2.6399996
     1 | γ = 0.69333315
     1 | p^2 = 0.7511107
     1 | s^2 = 0.7871538
     1 | u = -0.4932832
     2 | γ = -0.82671684
     2 | p^2 = 3.211054
     2 | s^2 = 0.20382473
     2 | u = 0.23642644
     3 | γ = 1.7502401
     3 | p^2 = 3.8475707
     3 | s^2 = 0.0
     3 | u = 0.0

     actual a
     0 | a = 6.6399994
     1 | a = 0.20004994
     2 | a = -0.5902904
     3 | a = 1.7502401
     4 | a = 0.0

     actual b
     0 | b^2 = 1.2704002, |b| = 1.1271203
     1 | b^2 = 3.1746697, |b| = 1.7817603
     2 | b^2 = 0.78423005, |b| = 0.88556767
     3 | b^2 = 0.0, |b| = 0.0

     iterating along diagonal
     0 | γ = 6.6399994
     0 | p^2 = 44.089592
     0 | s^2 = 0.028007062
     0 | u = 0.19156969
     1 | γ = 0.008480251
     1 | p^2 = 7.398681e-05
     1 | s^2 = 0.9999767
     1 | u = -0.58179665
     2 | γ = -0.008493781
     2 | p^2 = 3.095603
     2 | s^2 = 0.20212986
     2 | u = 0.35205892
     3 | γ = 1.3981812
     3 | p^2 = 2.4501615
     3 | s^2 = 0.0
     3 | u = 0.0

     actual a
     0 | a = 6.831569
     1 | a = -0.5733164
     2 | a = 0.34356514
     3 | a = 1.3981812
     4 | a = 0.0

     actual b
     0 | b^2 = 0.088915244, |b| = 0.2981866
     1 | b^2 = 3.8797426, |b| = 1.9697062
     2 | b^2 = 0.4952508, |b| = 0.7037406
     3 | b^2 = 0.0, |b| = 0.0

     iterating along diagonal
     0 | γ = 6.831569
     0 | p^2 = 46.670338
     0 | s^2 = 0.001901554
     0 | u = 0.011900405
     1 | γ = -0.5852168
     1 | p^2 = 0.3431312
     1 | s^2 = 0.9187446
     1 | u = -0.22201619
     2 | γ = 0.5655813
     2 | p^2 = 3.9367518
     2 | s^2 = 0.11174425
     2 | u = 0.21943916
     3 | γ = 1.178742
     3 | p^2 = 1.564226
     3 | s^2 = 0.0
     3 | u = 0.0

     actual a
     0 | a = 6.8434696
     1 | a = -0.807233
     2 | a = 0.7850205
     3 | a = 1.178742
     4 | a = 0.0

     actual b
     0 | b^2 = 0.0080300225, |b| = 0.08961039
     1 | b^2 = 4.0718784, |b| = 2.0178895
     2 | b^2 = 0.17479326, |b| = 0.41808283
     3 | b^2 = 0.0, |b| = 0.0

     iterating along diagonal
     0 | γ = 6.8434696
     0 | p^2 = 46.833076
     0 | s^2 = 0.00017143111
     0 | u = 0.0010347988
     1 | γ = -0.8082678
     1 | p^2 = 0.65340877
     1 | s^2 = 0.8617208
     1 | u = -0.020032683
     2 | γ = 0.8050532
     2 | p^2 = 4.686971
     2 | s^2 = 0.03595264
     2 | u = 0.07132267
     3 | γ = 1.1074194
     3 | p^2 = 1.2721134
     3 | s^2 = 0.0
     3 | u = 0.0

     actual a
     0 | a = 6.8445044
     1 | a = -0.8283005
     2 | a = 0.87637585
     3 | a = 1.1074194
     4 | a = 0.0

     actual b
     0 | b^2 = 0.00081006123, |b| = 0.028461576
     1 | b^2 = 4.1894836, |b| = 2.0468228
     2 | b^2 = 0.045735836, |b| = 0.21385938
     3 | b^2 = 0.0, |b| = 0.0

     iterating along diagonal
     0 | γ = 6.8445044
     0 | p^2 = 46.84724
     0 | s^2 = 1.7291248e-05
     0 | u = 0.00010402767
     1 | γ = -0.8284045
     1 | p^2 = 0.6862658
     1 | s^2 = 0.8592491
     1 | u = 0.041219354
     2 | γ = 0.8351565
     2 | p^2 = 4.955467
     2 | s^2 = 0.0091449665
     2 | u = 0.017764792
     3 | γ = 1.0896546
     3 | p^2 = 1.1983055
     3 | s^2 = 0.0
     3 | u = 0.0

     actual a
     0 | a = 6.8446083
     1 | a = -0.78718513
     2 | a = 0.8529213
     3 | a = 1.0896546
     4 | a = 0.0

     actual b
     0 | b^2 = 8.4307794e-05, |b| = 0.009181928
     1 | b^2 = 4.2972794, |b| = 2.072988
     2 | b^2 = 0.010958464, |b| = 0.104682684
     3 | b^2 = 0.0, |b| = 0.0

     */
    
    /*
     Before p^2 change:
     
     iterating along diagonal
     0 | γ = 2.0133333
     0 | p^2 = 4.053511
     0 | s^2 = 0.6894696
     0 | u = 2.3166175
     1 | γ = -0.9699513
     1 | p^2 = 3.0296726
     1 | s^2 = 0.4783129
     1 | u = -2.0455616
     2 | γ = -1.2611051
     2 | p^2 = 3.048544
     2 | s^2 = 0.21238251
     2 | u = -0.26783666
     3 | γ = 0.26783666
     3 | p^2 = 0.09108035
     3 | s^2 = 0.0
     3 | u = 0.0

     actual a
     0 | a = 6.3166175
     1 | a = -1.0288464
     2 | a = 0.4577248
     3 | a = 2.2545033
     4 | a = 0.0

     actual b
     0 | b^2 = 4.0040607, |b| = 2.001015
     1 | b^2 = 1.8513528, |b| = 1.3606442
     2 | b^2 = 0.019343873, |b| = 0.13908225
     3 | b^2 = 0.0, |b| = 0.0

     iterating along diagonal
     0 | γ = 4.0621142
     0 | p^2 = 16.500772
     0 | s^2 = 0.195274
     0 | u = 0.1520725
     1 | γ = -3.435422
     1 | p^2 = 14.666016
     1 | s^2 = 0.112085216
     1 | u = -0.58645236
     2 | γ = -1.2103261
     2 | p^2 = 1.6498083
     2 | s^2 = 0.011589042
     2 | u = -0.01402652
     3 | γ = 0.01402652
     3 | p^2 = 0.00019905006
     3 | s^2 = 0.0
     3 | u = 0.0

     actual a
     0 | a = 6.46869
     1 | a = -1.7673712
     2 | a = 1.0301507
     3 | a = 2.2685297
     4 | a = 0.0

     actual b
     0 | b^2 = 3.2254126, |b| = 1.7959434
     1 | b^2 = 0.18708728, |b| = 0.4325359
     2 | b^2 = 2.3067994e-06, |b| = 0.0015188152
     3 | b^2 = 0.0, |b| = 0.0

     iterating along diagonal
     0 | γ = 4.20016
     0 | p^2 = 17.641344
     0 | s^2 = 0.15457183
     0 | u = 0.025389807
     1 | γ = -4.0612907
     1 | p^2 = 19.509739
     1 | s^2 = 0.009498347
     1 | u = -0.0503381
     2 | γ = -1.1880409
     2 | p^2 = 1.424976
     2 | s^2 = 1.6188313e-06
     2 | u = -1.9232377e-06
     3 | γ = 1.9232377e-06
     3 | p^2 = 3.6988494e-12
     3 | s^2 = 0.0
     3 | u = 0.0

     actual a
     0 | a = 6.4940796
     1 | a = -1.8430994
     2 | a = 1.0804869
     3 | a = 2.2685316
     4 | a = 0.0

     actual b
     0 | b^2 = 3.0445745, |b| = 1.7448709
     1 | b^2 = 0.013534937, |b| = 0.11633975
     2 | b^2 = 5.9878132e-18, |b| = 2.4470008e-09
     3 | b^2 = 0.0, |b| = 0.0

     iterating along diagonal
     0 | γ = 4.225548
     0 | p^2 = 17.855253
     0 | s^2 = 0.14567462
     0 | u = 0.016594797
     1 | γ = -4.128226
     1 | p^2 = 19.948193
     1 | s^2 = 0.0006780444
     1 | u = -0.0036046673
     2 | γ = -1.18444
     2 | p^2 = 1.4038501
     2 | s^2 = 4.2652795e-18
     2 | u = -5.0519675e-18
     3 | γ = 5.0519675e-18
     3 | p^2 = 2.5522376e-35
     3 | s^2 = 0.0
     3 | u = 0.0

     actual a
     0 | a = 6.5106745
     1 | a = -1.8632991
     2 | a = 1.0840915
     3 | a = 2.2685316
     4 | a = 0.0

     actual b
     0 | b^2 = 2.907917, |b| = 1.7052616
     1 | b^2 = 0.0009518727, |b| = 0.030852433
     2 | b^2 = 0.0, |b| = 0.0
     3 | b^2 = 0.0, |b| = 0.0

     iterating along diagonal
     0 | γ = 4.2421427
     0 | p^2 = 17.995775
     0 | s^2 = 0.13911021
     0 | u = 0.015345523
     1 | γ = -4.1471763
     1 | p^2 = 19.97825
     1 | s^2 = 4.7643178e-05
     1 | u = -0.00025401515
     2 | γ = -1.184186
     2 | p^2 = 1.4023632
     2 | s^2 = 0.0
     2 | u = -0.0
     3 | γ = 0.0
     3 | p^2 = 0.0
     3 | s^2 = nan
     3 | u = nan

     actual a
     0 | a = 6.52602
     1 | a = -1.8788989
     2 | a = 1.0843456
     3 | a = nan
     4 | a = 0.0

     actual b
     0 | b^2 = 2.779311, |b| = 1.6671265
     1 | b^2 = 6.681304e-05, |b| = 0.008173924
     2 | b^2 = 0.0, |b| = 0.0
     3 | b^2 = 0.0, |b| = 0.0
     */
  }
}
