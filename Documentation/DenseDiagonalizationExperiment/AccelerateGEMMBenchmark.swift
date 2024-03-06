import XCTest
import Accelerate // Gate out once there are Swift kernels for CPUs without AMX.
import Numerics
import QuartzCore

// All matrices in this experiment are row-major. The individual wavefunctions
// are rows of the eigenvector matrix. However, LAPACK treats arguments as if
// they're column-major.
final class DenseDiagonalizationTests: XCTestCase {
  // MARK: - Algorithms
  
  // Multiplies two matrices in mixed precision, as a reference implementation
  // for faster methods.
  // - Requires row-major ordering.
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
  // - Requires row-major ordering.
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
  // - Row-major order: orthonormalizes the rows.
  // - Column-major order: orthonormalizes the columns.
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
  // - Row-major order: orthonormalizes the rows.
  // - Column-major order: orthonormalizes the columns.
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
  
  // Benchmark of Accelerate performance with skinny matrices.
  func testBlockSizePerformance() throws {
    // Tasks:
    // - Make a test to validate correctness of the Accelerate kernel.
    // - Benchmark matrices with height in { 512, 1024, 2048 }.
    // - Archive the functions above this, along with their updated comments.
    
    struct BenchmarkDescriptor {
      // The block size (n x n).
      var n: Int = 1
      
      // The serial duration of each trial (matrix height).
      var height: Int = 1
      
      // Whether to read matrix A transposed in memory.
      var transposeA: Bool = false
      
      // Whether to read matrix B transposed in memory.
      var transposeB: Bool = false
      
      // The number of samples to collect for the benchmark.
      var trialCount: Int = 3
    }
    
    var descriptor = BenchmarkDescriptor()
    descriptor.transposeA = true
    descriptor.transposeB = true
    descriptor.trialCount = 10
    
    // Record the data.
    let blockSizes: [Int] = [
      1, 2, 4, 8, 16, 32, 64, 128, 256, 512 , 1024
    ]
    var gflopsSummary: [Float] = []
    for blockSize in blockSizes {
      descriptor.n = blockSize
      print()
      
      var heights: [Int] = [
        4, 8, 16, 32, 64, 128, 256, 512, 1024 , 2048
      ]
      if blockSize == 1 {
        heights.removeFirst(4)
      }
      var gflopsMax: Float = .zero
      for height in heights {
        descriptor.height = height
        
        let latency = runBenchmark(descriptor: descriptor)
        let computeCost = descriptor.n * descriptor.n * descriptor.height
        let gflopsK = Float(computeCost) / Float(latency) / 1e9
        let gflopsRepr = String(format: "%.1f", gflopsK)
        let us = latency * 1e6
        let usRepr = String(format: "%.2f", us)
        print("n=\(descriptor.n)", terminator: ", ")
        print("h=\(descriptor.height)", terminator: ": ")
        print("\(gflopsRepr) GFLOPS/k", terminator: " ")
        print("(\(usRepr) μs)", terminator: "\n")
        
        gflopsMax = max(gflopsMax, gflopsK)
      }
      gflopsSummary.append(gflopsMax)
    }
    
    // Print the summary.
    print()
    print("summary")
    for blockSizeID in blockSizes.indices {
      let n = blockSizes[blockSizeID]
      let gflopsK = gflopsSummary[blockSizeID]
      let gflopsRepr = String(format: "%.1f", gflopsK)
      print("n=\(n)", terminator: ": ")
      print("\(gflopsRepr) GFLOPS/k", terminator: "\n")
    }
    
    // Run a benchmark with the specified configuration.
    // - returns: The minimum latency in seconds.
    func runBenchmark(descriptor: BenchmarkDescriptor) -> Float {
      let n = descriptor.n
      let height = descriptor.height
      let transposeA = descriptor.transposeA
      let transposeB = descriptor.transposeB
      let trialCount = descriptor.trialCount
      
      // These matrices are stored in row-major order.
      var matrixA = [Float](repeating: 0, count: height * n)
      var matrixB = [Float](repeating: 0, count: height * n)
      var expectedMatrixC = [Float](repeating: 0, count: n * n)
      for rowID in 0..<height {
        for columnID in 0..<n {
          let address = rowID * n + columnID
          let randomA = Float.random(in: -1...1)
          let randomB = Float.random(in: -1...1)
          matrixA[address] = randomA
          matrixB[address] = randomB
        }
      }
      for columnA in 0..<n {
        for columnB in 0..<n {
          var dotProduct: Double = .zero
          for rowID in 0..<height {
            let valueA = matrixA[rowID * n + columnA]
            let valueB = matrixB[rowID * n + columnB]
            dotProduct += Double(valueA) * Double(valueB)
          }
          
          // LAPACK returns C in reverse order (column-major). We should aim for
          // that format as well.
          let address = columnB * n + columnA
          expectedMatrixC[address] = Float(dotProduct)
        }
      }
      
      func transpose(matrix: [Float]) -> [Float] {
        var transposedMatrix = [Float](repeating: 0, count: n * height)
        for rowID in 0..<height {
          for columnID in 0..<n {
            let addressOriginal = rowID * n + columnID
            let addressTransposed = columnID * height + rowID
            transposedMatrix[addressTransposed] = matrix[addressOriginal]
          }
        }
        return transposedMatrix
      }
      var transposedMatrixA: [Float]
      var transposedMatrixB: [Float]
      if transposeA {
        transposedMatrixA = transpose(matrix: matrixA)
      } else {
        transposedMatrixA = matrixA
      }
      if transposeB {
        transposedMatrixB = matrixB
      } else {
        transposedMatrixB = transpose(matrix: matrixB)
      }
      
      var minLatency: Float = .greatestFiniteMagnitude
      var actualMatrixC = [Float](repeating: 0, count: n * n)
      for trialID in 0..<trialCount + 1 {
        var TRANSA = CChar(Character(transposeA ? "T" : "N").asciiValue!)
        var TRANSB = CChar(Character(transposeB ? "T" : "N").asciiValue!)
        var M: Int32 = Int32(n)
        var N: Int32 = Int32(n)
        var K: Int32 = Int32(height)
        var ALPHA: Float = 1
        var LDA: Int32 = transposeA ? Int32(height) : Int32(n)
        var BETA: Float = 0
        var LDB: Int32 = transposeB ? Int32(n) : Int32(height)
        var LDC: Int32 = Int32(n)
        transposedMatrixA.withContiguousStorageIfAvailable {
          let A = UnsafeMutablePointer(mutating: $0.baseAddress!)
          transposedMatrixB.withContiguousStorageIfAvailable {
            let B = UnsafeMutablePointer(mutating: $0.baseAddress!)
            actualMatrixC.withContiguousMutableStorageIfAvailable {
              let C = $0.baseAddress!
              if trialID < trialCount {
                let start = CACurrentMediaTime()
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
                let end = CACurrentMediaTime()
                let elapsedTime = Float(end - start)
                minLatency = min(minLatency, elapsedTime)
              } else {
                let start = CACurrentMediaTime()
                for _ in 0..<trialCount {
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
                let end = CACurrentMediaTime()
                let elapsedTime = Float(end - start) / Float(trialCount)
                minLatency = min(minLatency, elapsedTime)
              }
            }
          }
        }
        
        if trialID == 0 {
          var passRate: Float = .zero
          for columnA in 0..<n {
            for columnB in 0..<n {
              // LAPACK returns C in reverse order (column-major)
              let address = columnB * n + columnA
              let actual = actualMatrixC[address]
              let expected = expectedMatrixC[address]
              
              let error = 1e-5 * Float(height).squareRoot()
              if (actual - expected).magnitude < error {
                passRate += 1 / Float(n * n)
              }
              XCTAssertEqual(actual, expected, accuracy: error)
            }
          }
          XCTAssertGreaterThan(
            passRate, 0.99, "\n \(actualMatrixC) \n \(expectedMatrixC) \n")
        }
      }
      return minLatency
    }
  }
  
  // TODO: Create a custom direct diagonalization kernel. Start with the
  // simplest one known, then move on to more advanced solvers.
}
