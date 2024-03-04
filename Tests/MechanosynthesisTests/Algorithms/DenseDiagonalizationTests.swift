import XCTest
import Accelerate
import Numerics

// All matrices in this experiment are row-major. The individual wavefunctions
// are rows of the eigenvector matrix. However, LAPACK treats arguments as if
// they're column-major.

final class DenseDiagonalizationTests: XCTestCase {
  // MARK: - Algorithms
  
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
  
  // MARK: - Tests
  
  func testGramSchmidt() throws {
    // Generate a pre-determined matrix, test that the current version of
    // Gram-Schmidt normalizes them.
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
    // Generate a pre-determined matrix, test that the current version of
    // Gram-Schmidt normalizes them.
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
}
