import XCTest
import Accelerate // replace with a cross-platform import later
import Numerics
import QuartzCore // replace with a cross-platform import later

// All matrices in this experiment are row-major. The individual wavefunctions
// are rows of the eigenvector matrix. However, LAPACK treats arguments as if
// they're column-major.
final class DenseDiagonalizationTests: XCTestCase {
  // MARK: - Algorithms
  
  // Multiplies two matrices in mixed precision, as a reference implementation
  // for faster methods.
  static func matrixMultiply(
    matrixA: [Float], transposeA: Bool,
    matrixB: [Float], transposeB: Bool,
    matrixC: inout [Float], n: Int
  ) {
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
  }
  
  // Block version of matrix multiplication, utilizing the AMX coprocessor.
  static func blockMatrixMultiply(
    matrixA: [Float], transposeA: Bool,
    matrixB: [Float], transposeB: Bool,
    matrixC: inout [Float], n: Int
  ) {
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
  }
  
  // Orthonormalizes the matrix in mixed precision, as a reference
  // implementation for more approximate methods.
  static func gramSchmidtOrthonormalize(
    matrix: inout [Float], n: Int
  ) {
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
  }
  
  // Panel-factorized version of Gram-Schmidt. Its computation time scales
  // quadratically with the number of matrix elements (matrix K * panel count),
  // assuming infinite compute power.
  static func panelGramSchmidtOrthonormalize(
    matrix: inout [Float], n: Int, panelSize: Int
  ) {
    func normalize(electronID: Int) {
      var norm: Float = .zero
      for cellID in 0..<n {
        let value = matrix[electronID * n + cellID]
        norm += value * value
      }
      
      let normalizationFactor = 1 / norm.squareRoot()
      for cellID in 0..<n {
        var value = matrix[electronID * n + cellID]
        value *= normalizationFactor
        matrix[electronID * n + cellID] = value
      }
    }
    
    for electronID in 0..<n {
      // Normalize the vectors before taking dot products.
      normalize(electronID: electronID)
    }
    
    var dotProductMatrix = [Float](repeating: .zero, count: panelSize * n)
    var updateMatrix = [Float](repeating: .zero, count: panelSize * n)
    
    var panelStart = 0
    while panelStart < n {
      let panelEnd = min(panelStart + panelSize, n)
      
      // Determine the magnitude of components parallel to previous vectors.
      //  for electronID in panelStart..<panelEnd {
      //    let intraPanelID = electronID - panelStart
      //    for neighborID in 0..<panelStart {
      //      var dotProduct: Float = .zero
      //      for cellID in 0..<n {
      //        let value1 = matrix[electronID * n + cellID]
      //        let value2 = matrix[neighborID * n + cellID]
      //        dotProduct += value1 * value2
      //      }
      //      dotProductMatrix[intraPanelID * n + neighborID] = dotProduct
      //    }
      //  }
      do {
        var TRANSA = CChar(Character("T").asciiValue!)
        var TRANSB = CChar(Character("N").asciiValue!)
        var M: Int32 = Int32(panelStart)
        var N: Int32 = Int32(panelEnd - panelStart)
        var K: Int32 = Int32(n)
        var ALPHA: Float = 1
        var LDA: Int32 = Int32(n)
        var BETA: Float = 0
        var LDB: Int32 = Int32(n)
        var LDC: Int32 = Int32(n)
        matrix.withContiguousMutableStorageIfAvailable {
          let A = $0.baseAddress!
          let B = $0.baseAddress! + panelStart * n
          dotProductMatrix.withContiguousMutableStorageIfAvailable {
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
      
      // Negate all components parallel to previous vectors.
      //  for electronID in panelStart..<panelEnd {
      //    let intraPanelID = electronID - panelStart
      //    for cellID in 0..<n {
      //      var cellUpdate: Float = .zero
      //      for neighborID in panelStart..<electronID {
      //        let dotProduct = dotProductMatrix[intraPanelID * n + neighborID]
      //        let value2 = matrix[neighborID * n + cellID]
      //        cellUpdate -= dotProduct * value2
      //      }
      //      updateMatrix[intraPanelID * n + cellID] += cellUpdate
      //    }
      //  }
      do {
        var TRANSA = CChar(Character("N").asciiValue!)
        var TRANSB = CChar(Character("N").asciiValue!)
        var M: Int32 = Int32(n)
        var N: Int32 = Int32(panelEnd - panelStart)
        var K: Int32 = Int32(panelStart)
        var ALPHA: Float = -1
        var LDA: Int32 = Int32(n)
        var BETA: Float = 0
        var LDB: Int32 = Int32(n)
        var LDC: Int32 = Int32(n)
        matrix.withContiguousMutableStorageIfAvailable {
          let A = $0.baseAddress!
          dotProductMatrix.withContiguousMutableStorageIfAvailable {
            let B = $0.baseAddress!
            updateMatrix.withContiguousMutableStorageIfAvailable {
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
      }
      
      for electronID in panelStart..<panelEnd {
        let intraPanelID = electronID - panelStart
        
        // Determine the magnitude of components parallel to previous vectors.
        for neighborID in panelStart..<electronID {
          var dotProduct: Float = .zero
          for cellID in 0..<n {
            let value1 = matrix[electronID * n + cellID]
            let value2 = matrix[neighborID * n + cellID]
            dotProduct += value1 * value2
          }
          dotProductMatrix[intraPanelID * n + neighborID] = dotProduct
        }
        
        // Negate all components parallel to previous vectors.
        for cellID in 0..<n {
          var cellUpdate: Float = .zero
          for neighborID in panelStart..<electronID {
            let dotProduct = dotProductMatrix[intraPanelID * n + neighborID]
            let value2 = matrix[neighborID * n + cellID]
            cellUpdate -= dotProduct * value2
          }
          updateMatrix[intraPanelID * n + cellID] += cellUpdate
        }
        
        // Add the update to the vector.
        for cellID in 0..<n {
          var value1 = matrix[electronID * n + cellID]
          let cellUpdate = updateMatrix[intraPanelID * n + cellID]
          value1 += cellUpdate
          matrix[electronID * n + cellID] = value1
        }
        
        // Rescale the remaining components to unit length.
        normalize(electronID: electronID)
      }
      panelStart += panelSize
    }
  }
  
  // Direct diagonalization using divide-and-conquer from LAPACK.
  // - eigenvalues: n-element array
  // - eigenvalues: n x n matrix
  static func diagonalize(
    matrix: [Float], n: Int
  ) -> (eigenvalues: [Float], eigenvectors: [Float]) {
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
    
    return (W, A)
  }
  
  // Direct derivation of eigenvalues using 2-stage tridiagonalization.
  static func findEigenvalues(
    matrix: [Float], n: Int
  ) -> [Float] {
    var JOBZ = CChar(Character("N").asciiValue!)
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
    
    return W
  }
  
  // MARK: - Tests
  
  func testMatrixMultiply() throws {
    let a: Float = 1.5
    let b: Float = 2.5
    let c: Float = 3.5
    let d: Float = 4.5
    let ac = a * c
    let ad = a * d
    let bc = b * c
    let bd = b * d
    let matrixA: [Float] = [
      a, b, 0, 0,
      0, a, b, 0,
      0, 0, a, b,
      0, 0, 0, a,
    ]
    let matrixB: [Float] = [
      c, d, 0, 0,
      0, c, d, 0,
      0, 0, c, d,
      0, 0, 0, c,
    ]
    var matrixC = [Float](repeating: 0, count: 4 * 4)
    
    typealias MultiplyFunction = (
      [Float], Bool,
      [Float], Bool,
      inout [Float], Int
    ) -> Void
    var functions: [MultiplyFunction] = []
    functions.append(Self.matrixMultiply(
      matrixA:transposeA:matrixB:transposeB:matrixC:n:))
    functions.append(Self.blockMatrixMultiply(
      matrixA:transposeA:matrixB:transposeB:matrixC:n:))
    
    for multiplyFunction in functions {
      matrixC = [Float](repeating: 0, count: 4 * 4)
      multiplyFunction(
        matrixA, false,
        matrixB, false,
        &matrixC, 4)
      XCTAssertEqual(matrixC, [
        ac, ad + bc, bd, 0,
        0, ac, ad + bc, bd,
        0, 0, ac, ad + bc,
        0, 0, 0, ac,
      ])
      
      matrixC = [Float](repeating: 0, count: 4 * 4)
      multiplyFunction(
        matrixA, false,
        matrixB, true,
        &matrixC, 4)
      XCTAssertEqual(matrixC, [
        ac + bd, bc, 0, 0,
        ad, ac + bd, bc, 0,
        0, ad, ac + bd, bc,
        0, 0, ad, ac,
      ])
      
      matrixC = [Float](repeating: 0, count: 4 * 4)
      multiplyFunction(
        matrixA, true,
        matrixB, false,
        &matrixC, 4)
      XCTAssertEqual(matrixC, [
        ac, ad, 0, 0,
        bc, bd + ac, ad, 0,
        0, bc, bd + ac, ad,
        0, 0, bc, bd + ac,
      ])
      
      matrixC = [Float](repeating: 0, count: 4 * 4)
      multiplyFunction(
        matrixA, true,
        matrixB, true,
        &matrixC, 4)
      XCTAssertEqual(matrixC, [
        ac, 0, 0, 0,
        ad + bc, ac, 0, 0,
        bd, ad + bc, ac, 0,
        0, bd, ad + bc, ac,
      ])
    }
  }
  
  func testGramSchmidt() throws {
    var matrix: [Float] = [
      7, 6, 5, 4, 3, 2, 1,
      6, 7, 5, 4, 3, 2, 1,
      2, -1, -1, 1, 2, 6, 8,
      0.1, 0.2, 0.5, 0.5, 0.1, 0.2, 0.5,
      -0.1, 0.2, -0.5, 0.5, -0.1, 0.2, -0.5,
      -1, -2, -3, -5, -7, -9, -10,
      69, 23, 9, -48, 7, 1, 9,
    ]
    guard matrix.count == 49 else {
      fatalError("Not a 7x7 matrix.")
    }
    Self.gramSchmidtOrthonormalize(matrix: &matrix, n: 7)
    
    var overlapMatrix = [Float](repeating: .zero, count: 49)
    for electronID in 0..<7 {
      for neighborID in 0..<7 {
        var overlap: Float = .zero
        for cellID in 0..<7 {
          let value1 = matrix[electronID * 7 + cellID]
          let value2 = matrix[neighborID * 7 + cellID]
          overlap += value1 * value2
        }
        overlapMatrix[electronID * 7 + neighborID] = overlap
      }
    }
    
    for electronID in 0..<7 {
      for neighborID in 0..<7 {
        let value = overlapMatrix[electronID * 7 + neighborID]
        if electronID == neighborID {
          XCTAssertEqual(value, 1.0, accuracy: 1e-5)
        } else {
          XCTAssertEqual(value, 0.0, accuracy: 1e-5)
        }
      }
    }
  }
  
  func testPanelGramSchmidt() throws {
    var matrix: [Float] = [
      7, 6, 5, 4, 3, 2, 1,
      6, 7, 5, 4, 3, 2, 1,
      2, -1, -1, 1, 2, 6, 8,
      0.1, 0.2, 0.5, 0.5, 0.1, 0.2, 0.5,
      -0.1, 0.2, -0.5, 0.5, -0.1, 0.2, -0.5,
      -1, -2, -3, -5, -7, -9, -10,
      69, 23, 9, -48, 7, 1, 9,
    ]
    guard matrix.count == 49 else {
      fatalError("Not a 7x7 matrix.")
    }
    Self.panelGramSchmidtOrthonormalize(matrix: &matrix, n: 7, panelSize: 4)
    
    var overlapMatrix = [Float](repeating: .zero, count: 49)
    for electronID in 0..<7 {
      for neighborID in 0..<7 {
        var overlap: Float = .zero
        for cellID in 0..<7 {
          let value1 = matrix[electronID * 7 + cellID]
          let value2 = matrix[neighborID * 7 + cellID]
          overlap += value1 * value2
        }
        overlapMatrix[electronID * 7 + neighborID] = overlap
      }
    }
    
    for electronID in 0..<7 {
      for neighborID in 0..<7 {
        let value = overlapMatrix[electronID * 7 + neighborID]
        if electronID == neighborID {
          XCTAssertEqual(value, 1.0, accuracy: 1e-4)
        } else {
          XCTAssertEqual(value, 0.0, accuracy: 1e-4)
        }
      }
    }
  }
  
  // Test the performance of eigensolvers in LAPACK.
  func testDiagonalization() throws {
    let checkingEigenvalueCorrectness = true
    
    func benchmarkProblemSize(n: Int, trialCount: Int) {
      func createEigenvalue(electronID: Int) -> Float {
        let position = Float(n / 2 - electronID)
        var sign: Float
        if position > 0 {
          sign = 1
        } else if position == 0 {
          sign = 0
        } else {
          sign = -1
        }
        
        let powerPart = Float.pow(1.02, position.magnitude)
        return sign * powerPart
      }
      
      var Λ = [Float](repeating: 0, count: n * n)
      var ΣT = [Float](repeating: 0, count: n * n)
      for electronID in 0..<n {
        let eigenvalue = createEigenvalue(electronID: electronID)
        let address = electronID * n + electronID
        Λ[address] = eigenvalue
      }
      for entry in ΣT.indices {
        ΣT[entry] = Float.random(in: -1...1)
      }
      Self.gramSchmidtOrthonormalize(matrix: &ΣT, n: n)
      
      for electronID in 0..<n {
        var norm: Float = .zero
        for cellID in 0..<n {
          let address = electronID * n + cellID
          norm += ΣT[address] * ΣT[address]
        }
        XCTAssertEqual(norm, 1, accuracy: 1e-3)
      }
      
      var ΛΣT = [Float](repeating: 0, count: n * n)
      var ΣΛΣT = [Float](repeating: 0, count: n * n)
      Self.blockMatrixMultiply(
        matrixA: Λ, transposeA: false,
        matrixB: ΣT, transposeB: false,
        matrixC: &ΛΣT, n: n)
      Self.blockMatrixMultiply(
        matrixA: ΣT, transposeA: true,
        matrixB: ΛΣT, transposeB: false,
        matrixC: &ΣΛΣT, n: n)
      
      let A = ΣΛΣT
      var hamiltonianPsi = [Float](repeating: 0, count: n * n)
      Self.blockMatrixMultiply(
        matrixA: ΣT, transposeA: false,
        matrixB: A, transposeB: true,
        matrixC: &hamiltonianPsi, n: n)
      
      do {
        var eigenvaluePassRate: Float = .zero
        for electronID in 0..<n {
          var λ: Float = .zero
          for cellID in 0..<n {
            let address = electronID * n + cellID
            λ += ΣT[address] * hamiltonianPsi[address]
          }
          
          let address = electronID * n + electronID
          let eigenvalue = Λ[address]
          let error = max(
            1, createEigenvalue(electronID: electronID).magnitude)
          
          if (eigenvalue - λ).magnitude < 1e-3 * error {
           eigenvaluePassRate += 1 / Float(n)
          }
        }
        XCTAssertGreaterThan(eigenvaluePassRate, 0.90)
      }
      
      var minEigenvalueTime: Float = .greatestFiniteMagnitude
      var minEigenvectorTime: Float = .greatestFiniteMagnitude
      for trialID in 0..<trialCount {
        do {
          let start = CACurrentMediaTime()
          let eigenvalues = Self.findEigenvalues(matrix: A, n: n)
          let end = CACurrentMediaTime()
          minEigenvalueTime = min(minEigenvalueTime, Float(end - start))
          
          if trialID == 0 {
            // Validate correctness of eigenvalues.
            var eigenvaluePassRate: Float = .zero
            for electronID in 0..<n {
              let expected = Λ[electronID * n + electronID]
              let actual = eigenvalues[n - 1 - electronID]
              let error = max(
                1, createEigenvalue(electronID: electronID).magnitude)
              
              if (actual - expected).magnitude < 1e-3 * error {
                eigenvaluePassRate += 1 / Float(n)
              }
            }
            if checkingEigenvalueCorrectness {
              XCTAssertGreaterThan(eigenvaluePassRate, 0.90)
            }
          }
        }
        
        do {
          let start = CACurrentMediaTime()
          let (eigenvalues, eigenvectors) = Self.diagonalize(matrix: A, n: n)
          let end = CACurrentMediaTime()
          minEigenvectorTime = min(minEigenvectorTime, Float(end - start))
          
          if trialID == 0 {
            // Validate correctness of eigenvalues.
            var eigenvaluePassRate: Float = .zero
            for electronID in 0..<n {
              let expected = Λ[electronID * n + electronID]
              let actual = eigenvalues[n - 1 - electronID]
              let error = max(
                1, createEigenvalue(electronID: electronID).magnitude)
              
              if (actual - expected).magnitude < 1e-3 * error {
                eigenvaluePassRate += 1 / Float(n)
              }
            }
            XCTAssertGreaterThan(eigenvaluePassRate, 0.90)
            
            // Validate correctness of eigenvectors.
            var eigenvectorPassRate: Float = .zero
            for electronID in 0..<n {
              var dotProduct: Float = .zero
              for cellID in 0..<n {
                let address1 = electronID * n + cellID
                let address2 = (n - 1 - electronID) * n + cellID
                let expected = ΣT[address1]
                let actual = eigenvectors[address2]
                dotProduct += expected * actual
              }
              if (dotProduct.magnitude - 1).magnitude < 1e-3 {
                eigenvectorPassRate += 1 / Float(n)
              }
            }
            XCTAssertGreaterThan(eigenvectorPassRate, 0.90)
          }
        }
      }
      
      #if false
      minEigenvectorTime -= minEigenvalueTime
      let minOverallTime = minEigenvalueTime + minEigenvectorTime
      
      let latencyΛ = Int(minEigenvalueTime * 1e6)
      let latencyΣ = Int(minEigenvectorTime * 1e6)
      let latencyΛΣ = Int(minOverallTime * 1e6)
      func createGFLOPSk(time: Float) -> String {
        let computeCost = Float(n * n * n)
        let speed = computeCost / time / 1e9
        return String(format: "%.2f", speed)
      }
      let gflopsΛ = createGFLOPSk(time: minEigenvalueTime)
      let gflopsΣ = createGFLOPSk(time: minEigenvectorTime)
      let gflopsΛΣ = createGFLOPSk(time: minOverallTime)
      print("n=\(n)", terminator: " ")
      print("| latency (μs)", terminator: " ")
      print("Λ=\(latencyΛ) Σ=\(latencyΣ) ΛΣ=\(latencyΛΣ)", terminator: " ")
      print("| GFLOPS/k", terminator: " ")
      print("Λ=\(gflopsΛ) Σ=\(gflopsΣ) ΛΣ=\(gflopsΛΣ)", terminator: " ")
      print()
      #endif
    }
    
    for n in 1...20 {
      benchmarkProblemSize(n: n, trialCount: 3)
    }
    #if false
    for nTenth in 3...10 {
      benchmarkProblemSize(n: 10 * nTenth, trialCount: 3)
    }
    benchmarkProblemSize(n: 125, trialCount: 3)
    benchmarkProblemSize(n: 200, trialCount: 3)
    benchmarkProblemSize(n: 250, trialCount: 3)
    benchmarkProblemSize(n: 300, trialCount: 3)
    benchmarkProblemSize(n: 400, trialCount: 3)
    benchmarkProblemSize(n: 500, trialCount: 3)
    benchmarkProblemSize(n: 600, trialCount: 3)
    benchmarkProblemSize(n: 750, trialCount: 3)
    benchmarkProblemSize(n: 800, trialCount: 3)
    benchmarkProblemSize(n: 1000, trialCount: 3)
    #endif
  }
}
