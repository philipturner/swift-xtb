import XCTest
import Mechanosynthesis
import Numerics

import Accelerate

final class LinearSolverTests: XCTestCase {
  // Creates a 1D Laplacian operator with 2nd order accuracy.
  static func laplacian(h: Float) -> ([Float]) -> [Float] {
    return { potential in
      var chargeDensity: [Float] = []
      
      for cellID in potential.indices {
        var left: Float = .zero
        var center: Float = .zero
        var right: Float = .zero
        
        if potential.indices.contains(cellID - 1) {
          left = potential[cellID - 1]
        }
        if potential.indices.contains(cellID) {
          center = potential[cellID]
        }
        if potential.indices.contains(cellID + 1) {
          right = potential[cellID + 1]
        }
        
        let leftDerivative = (center - left) / h
        let rightDerivative = (right - center) / h
        let doubleDerivative = (rightDerivative - leftDerivative) / h
        chargeDensity.append(doubleDerivative)
      }
      return chargeDensity
    }
  }
  
  // LAPACK solver as a reference implementation.
  func testDirectMatrixMethod() throws {
    // The input matrix must be column-major. Code using this LAPACK function
    // will likely need a second function that transposes the matrix.
    var coefficients: [Float] = []
    coefficients += [Float(1), Float(1), Float(1)]
    coefficients += [Float(4), Float(2), Float(1)]
    coefficients += [Float(9), Float(3), Float(1)]
    
    var rightHandSide: [Float] = []
    rightHandSide += [Float(9), Float(3), Float(1)]
    
    var N: Int32 = 3
    var NRHS: Int32 = 1
    var A: [Float] = coefficients
    var LDA: Int32 = 3
    var IPIV: [Int32] = .init(repeating: .zero, count: 3)
    var B: [Float] = rightHandSide
    var LDB: Int32 = 3
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
      }
    }
    XCTAssertEqual(INFO, 0, "Linear solver failed.")
    
    let X = B
    XCTAssertEqual(X[0], 0.000, accuracy: 1e-3)
    XCTAssertEqual(X[1], 0.000, accuracy: 1e-3)
    XCTAssertEqual(X[2], 1.000, accuracy: 1e-3)
  }
  
  func testSteepestDescent() throws {
    
  }
}
