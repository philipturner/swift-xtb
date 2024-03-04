import XCTest
import Accelerate
import Numerics

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
  
  // Custom eigensolver using iterative techniques.
  // - eigenvalues: n-element array
  // - eigenvalues: n x n matrix
  static func customDiagonalize(
    matrix H: [Float], n: Int
  ) -> (eigenvalues: [Float], eigenvectors: [Float]) {
    // Initialize the wavefunctions to the identity matrix.
    var E = [Float](repeating: .zero, count: n)
    var Ψ = [Float](repeating: .zero, count: n * n)
    for electronID in 0..<n {
      let address = electronID * n + electronID
      Ψ[address] = 1
    }
    
    func updateEnergies(HΨ: [Float]) {
      for electronID in 0..<n {
        var ΨHΨ: Float = .zero
        for cellID in 0..<n {
          let address = electronID * n + cellID
          ΨHΨ += Ψ[address] * HΨ[address]
        }
        E[electronID] = ΨHΨ
        print(ΨHΨ, terminator: ", ")
      }
      print()
    }
    
    func sortEigenpairs(rule: (Float, Float) -> Bool) {
      var sortedE: [SIMD2<Float>] = []
      for electronID in 0..<n {
        let key = E[electronID]
        let value = Float(electronID)
        sortedE.append(SIMD2(key, value))
      }
      sortedE.sort(by: { rule($0.x, $1.x) })
      
      var newE = Array(repeating: Float.zero, count: E.count)
      var newΨ = Array(repeating: Float.zero, count: Ψ.count)
      for newElectronID in 0..<n {
        let oldElectronID = Int(sortedE[newElectronID][1])
        newE[newElectronID] = E[oldElectronID]
        
        for cellID in 0..<n {
          let oldAddress = oldElectronID * n + cellID
          let newAddress = newElectronID * n + cellID
          newΨ[newAddress] = Ψ[oldAddress]
        }
      }
      E = newE
      Ψ = newΨ
    }
    
    print("diagonalization")
    for iterationID in 0..<60 {
      print("iteration \(iterationID)", terminator: ": ")
      
      // MARK: - Power Iteration (Reference Implementation)
      #if false
      // Multiply by the Hamiltonian and query the eigenvalues.
      var HΨ = [Float](repeating: .zero, count: n * n)
      Self.matrixMultiply(
        matrixA: Ψ, transposeA: false,
        matrixB: H, transposeB: true,
        matrixC: &HΨ, n: n)
      updateEnergies(HΨ: HΨ)
      
      // Refine each eigenvector in isolation.
      for electronID in 0..<n {
        for cellID in 0..<n {
          let address = electronID * n + cellID
          Ψ[address] = HΨ[address]
        }
      }
      
      // Sort the eigenpairs for the next iteration.
      sortEigenpairs(rule: { $0.magnitude > $1.magnitude })
      
      // Form an orthonormal set.
      Self.gramSchmidtOrthonormalize(matrix: &Ψ, n: n)
      
      // MARK: - Steepest Descent
      #else
      // Multiply by the Hamiltonian and query the eigenvalues.
      var HΨ = [Float](repeating: .zero, count: n * n)
      Self.blockMatrixMultiply(
        matrixA: Ψ, transposeA: false,
        matrixB: H, transposeB: true,
        matrixC: &HΨ, n: n)
      updateEnergies(HΨ: HΨ)
      
      // Refine each eigenvector in isolation.
      for _ in 0..<5 {
        var E = [Float](repeating: 0, count: n)
        var norm = [Float](repeating: 0, count: n)
        var r = [Float](repeating: 0, count: n * n)
        for electronID in 0..<n {
          var currentE: Float = .zero
          var currentNorm: Float = .zero
          for cellID in 0..<n {
            let address = electronID * n + cellID
            let valueΨ = Ψ[address]
            let valueHΨ = HΨ[address]
            currentE += valueΨ * valueHΨ
            currentNorm += valueΨ * valueΨ
          }
          E[electronID] = currentE
          norm[electronID] = currentNorm
          
          let rayleigh = currentE / currentNorm
          for cellID in 0..<n {
            let address = electronID * n + cellID
            let valueΨ = Ψ[address]
            let valueHΨ = HΨ[address]
            r[address] = valueHΨ - rayleigh * valueΨ
          }
        }
        
        // Multiply the residual by the Hamiltonian.
        var Hr = [Float](repeating: .zero, count: n * n)
        Self.blockMatrixMultiply(
          matrixA: r, transposeA: false,
          matrixB: H, transposeB: true,
          matrixC: &Hr, n: n)
        
        // Construct the timesteps.
        var λ = [Float](repeating: .zero, count: n)
        for electronID in 0..<n {
          var rr: Float = .zero
          var Ψr: Float = .zero
          var rHr: Float = .zero
          var ΨHr: Float = .zero
          for cellID in 0..<n {
            let address = electronID * n + cellID
            rr += r[address] * r[address]
            Ψr += Ψ[address] * r[address]
            rHr += r[address] * Hr[address]
            ΨHr += Ψ[address] * Hr[address]
          }
          
          let (m0, m1, m2, m3) = (rr, Ψr, rHr, ΨHr)
          let (m4, m5) = (E[electronID], norm[electronID])
          let ca = (m0 * m3 - m2 * m1)
          let cb = (m5 * m2 - m4 * m0)
          let cc = (m4 * m1 - m3 * m5)
          let determinant = cb + (cb * cb - 4 * ca * cc).squareRoot()
          if determinant.magnitude > 1e-15 {
            λ[electronID] = 2 * cc / determinant
          }
          
          if rr.squareRoot() < 1e-5 {
            λ[electronID] = 0
          }
        }
        for electronID in 0..<n {
          let timeStep = λ[electronID]
          for cellID in 0..<n {
            let address = electronID * n + cellID
            Ψ[address] += timeStep * r[address]
            HΨ[address] += timeStep * Hr[address]
          }
        }
      }
      
      // If a contiguous set of the lowest 'm' eigenvectors have all
      // converged, we can omit them from the O(n^3) orthogonalization
      // and the expensive SD iterations. Simply mask out the lower
      // portion of the matrix, making it a dense linear algebra kernel.
      sortEigenpairs(rule: { $0 < $1 })
      
      // Form an orthonormal set.
      Self.panelGramSchmidtOrthonormalize(
        matrix: &Ψ, n: n, panelSize: 4)
      #endif
    }
    
    sortEigenpairs(rule: { $0 > $1 })
    
    return (E, Ψ)
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
  
  // Test the eigensolvers with regular eigenvalues.
  func testEigendecomposition() throws {
    // Ψ is already stored in a transposed format.
    var Ψ: [Float] = [
      7, 6, 5, 4, 3, 2, 1,
      6, 7, 5, 4, 3, 2, 1,
      2, -1, -1, 1, 2, 6, 8,
      0.1, 0.2, 0.5, 0.5, 0.1, 0.2, 0.5,
      -0.1, 0.2, -0.5, 0.5, -0.1, 0.2, -0.5,
      -1, -2, -3, -5, -7, -9, -10,
      69, 23, 9, -48, 7, 1, 9,
    ]
    Self.gramSchmidtOrthonormalize(matrix: &Ψ, n: 7)
    
    print("eigenvectors")
    for electronID in 0..<7 {
      for cellID in 0..<7 {
        let value = Ψ[electronID * 7 + cellID]
        print(value, terminator: ", ")
      }
      print()
    }
    
    // Alternative eigenspectrum for testing the power method:
    // 4, 3, 2, 1, 0.1, -1.3, -2.3
    let eigenvalues: [Float] = [
      4, 3, 2, 1, 0, -1, -2
    ]
    var Λ = [Float](repeating: 0, count: 7 * 7)
    for i in 0..<7 {
      let address = i * 7 + i
      let value = eigenvalues[i]
      Λ[address] = value
    }
    
    print("eigenvalues")
    for electronID in 0..<7 {
      for neighborID in 0..<7 {
        let value = Λ[electronID * 7 + neighborID]
        print(value, terminator: ", ")
      }
      print()
    }
    
    var ΛΣT = [Float](repeating: 0, count: 7 * 7)
    Self.matrixMultiply(
      matrixA: Λ, transposeA: false,
      matrixB: Ψ, transposeB: false,
      matrixC: &ΛΣT, n: 7)
    var H = [Float](repeating: 0, count: 7 * 7)
    Self.matrixMultiply(
      matrixA: Ψ, transposeA: true,
      matrixB: ΛΣT, transposeB: false,
      matrixC: &H, n: 7)
    
    print("hamiltonian")
    for electronID in 0..<7 {
      for neighborID in 0..<7 {
        let value = H[electronID * 7 + neighborID]
        print(value, terminator: ", ")
      }
      print()
    }
    
    var HΨ = [Float](repeating: 0, count: 7 * 7)
    Self.matrixMultiply(
      matrixA: Ψ, transposeA: false,
      matrixB: H, transposeB: true,
      matrixC: &HΨ, n: 7)
    
    for electronID in 0..<7 {
      var rNorm: Float = .zero
      let expectedEigenvalue = eigenvalues[electronID]
      for cellID in 0..<7 {
        var value = HΨ[electronID * 7 + cellID]
        value -= expectedEigenvalue * Ψ[electronID * 7 + cellID]
        rNorm += value * value
      }
      rNorm.formSquareRoot()
      XCTAssertLessThan(rNorm, 1e-4)
    }
    
    let (lapackEigenvalues, lapackEigenvectors) = Self.diagonalize(
      matrix: H, n: 7)
    print("lapack eigenvalues")
    for electronID in 0..<7 {
      let expected = eigenvalues[6 - electronID]
      let actual = lapackEigenvalues[electronID]
      XCTAssertEqual(actual, expected, accuracy: 1e-5)
      print(actual, terminator: ", ")
    }
    print()
    
    print("lapack eigenvectors")
    for electronID in 0..<7 {
      var dotProduct: Float = .zero
      for cellID in 0..<7 {
        let reversedID = 6 - electronID
        let actual = lapackEigenvectors[reversedID * 7 + cellID]
        let expected = Ψ[electronID * 7 + cellID]
        print(actual, terminator: ", ")
        dotProduct += actual * expected
      }
      XCTAssertEqual(dotProduct.magnitude, 1, accuracy: 1e-5)
      print("[\(dotProduct)]")
    }
    
    let (customEigenvalues, customEigenvectors) = Self.customDiagonalize(
      matrix: H, n: 7)
    print("custom eigenvalues")
    for electronID in 0..<7 {
      let expected = eigenvalues[electronID]
      let actual = customEigenvalues[electronID]
      print(actual, terminator: ", ")
    }
    print()
    
    print("custom eigenvectors")
    for electronID in 0..<7 {
      var dotProduct: Float = .zero
      for cellID in 0..<7 {
        let actual = customEigenvectors[electronID * 7 + cellID]
        let expected = Ψ[electronID * 7 + cellID]
        print(actual, terminator: ", ")
        dotProduct += actual * expected
      }
      print("[\(dotProduct)]")
    }
  }
  
  // Test the eigensolvers with degenerate eigenvalue clusters.
  func testDegenerateEigenspaces() throws {
    // Ψ is already stored in a transposed format.
    var Ψ: [Float] = [
      7, 6, 5, 4, 3, 2, 1,
      6, 7, 5, 4, 3, 2, 1,
      2, -1, -1, 1, 2, 6, 8,
      0.1, 0.2, 0.5, 0.5, 0.1, 0.2, 0.5,
      -0.1, 0.2, -0.5, 0.5, -0.1, 0.2, -0.5,
      -1, -2, -3, -5, -7, -9, -10,
      69, 23, 9, -48, 7, 1, 9,
    ]
    Self.gramSchmidtOrthonormalize(matrix: &Ψ, n: 7)
    
    print("eigenvectors")
    for electronID in 0..<7 {
      for cellID in 0..<7 {
        let value = Ψ[electronID * 7 + cellID]
        print(value, terminator: ", ")
      }
      print()
    }
    
    // Cluster 1: close, but not equal eigenvalues
    // Cluster 2: pathological zero eigenvalue
    // Cluster 3: equal eigenvalues, meaning infinitely many eigenvectors
    let eigenvalues: [Float] = [
      4, 3.01, 3, 2.99, 0, -2, -2
    ]
    var Λ = [Float](repeating: 0, count: 7 * 7)
    for i in 0..<7 {
      let address = i * 7 + i
      let value = eigenvalues[i]
      Λ[address] = value
    }
    
    print("eigenvalues")
    for electronID in 0..<7 {
      for neighborID in 0..<7 {
        let value = Λ[electronID * 7 + neighborID]
        print(value, terminator: ", ")
      }
      print()
    }
    
    var ΛΣT = [Float](repeating: 0, count: 7 * 7)
    Self.matrixMultiply(
      matrixA: Λ, transposeA: false,
      matrixB: Ψ, transposeB: false,
      matrixC: &ΛΣT, n: 7)
    var H = [Float](repeating: 0, count: 7 * 7)
    Self.matrixMultiply(
      matrixA: Ψ, transposeA: true,
      matrixB: ΛΣT, transposeB: false,
      matrixC: &H, n: 7)
    
    print("hamiltonian")
    for electronID in 0..<7 {
      for neighborID in 0..<7 {
        let value = H[electronID * 7 + neighborID]
        print(value, terminator: ", ")
      }
      print()
    }
    
    var HΨ = [Float](repeating: 0, count: 7 * 7)
    Self.matrixMultiply(
      matrixA: Ψ, transposeA: false,
      matrixB: H, transposeB: true,
      matrixC: &HΨ, n: 7)
    
    for electronID in 0..<7 {
      var rNorm: Float = .zero
      let expectedEigenvalue = eigenvalues[electronID]
      for cellID in 0..<7 {
        var value = HΨ[electronID * 7 + cellID]
        value -= expectedEigenvalue * Ψ[electronID * 7 + cellID]
        rNorm += value * value
      }
      rNorm.formSquareRoot()
      XCTAssertLessThan(rNorm, 1e-4)
    }
    
    let (lapackEigenvalues, lapackEigenvectors) = Self.diagonalize(
      matrix: H, n: 7)
    print("lapack eigenvalues")
    for electronID in 0..<7 {
      let expected = eigenvalues[6 - electronID]
      let actual = lapackEigenvalues[electronID]
      XCTAssertEqual(actual, expected, accuracy: 1e-4)
      print(actual, terminator: ", ")
    }
    print()
    print()
    
    print("lapack eigenvectors")
    for electronID in 0..<7 {
      var dotProduct: Float = .zero
      for cellID in 0..<7 {
        let reversedID = 6 - electronID
        let actual = lapackEigenvectors[reversedID * 7 + cellID]
        let expected = Ψ[electronID * 7 + cellID]
        print(actual, terminator: ", ")
        dotProduct += actual * expected
      }
      if electronID != 5, electronID != 6 {
        XCTAssertEqual(dotProduct.magnitude, 1, accuracy: 1e-5)
      }
      print("[\(dotProduct)]")
    }
  }
}
