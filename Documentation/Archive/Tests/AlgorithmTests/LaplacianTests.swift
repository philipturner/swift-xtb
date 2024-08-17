import XCTest
import LinearAlgebra

final class LaplacianTests: XCTestCase {
  // Directly invert the Laplacian for a 4-cell domain.
  func testLaplacianInverse4x4() throws {
    // The matrices are in row-major order.
    func multiply4x4(_ A: [Float], _ B: [Float]) -> [Float] {
      LinearAlgebraUtilities.matrixMultiply(matrixA: A, matrixB: B, n: 4)
    }
    func transpose4x4(_ A: [Float]) -> [Float] {
      var output = [Float](repeating: 0, count: 16)
      for m in 0..<4 {
        for n in 0..<4 {
          let value = A[n * 4 + m]
          output[m * 4 + n] = value
        }
      }
      return output
    }
    
    // Create a 4x4 matrix that represents an aperiodic Laplacian.
    var laplacian: [Float] = []
    laplacian += [-2, 1, 0, 0]
    laplacian += [1, -2, 1, 0]
    laplacian += [0, 1, -2, 1]
    laplacian += [0, 0, 1, -2]
    
    // Hard-code the results of a 4x4 inversion with Accelerate.
    var laplacianInverse: [Float] = []
    laplacianInverse += [-0.8, -0.6, -0.4, -0.2]
    laplacianInverse += [-0.6, -1.2, -0.8, -0.4]
    laplacianInverse += [-0.4, -0.8, -1.2, -0.6]
    laplacianInverse += [-0.2, -0.4, -0.6, -0.8]
    
    // Evaluate the product of these two matrices.
    let product = multiply4x4(laplacian, laplacianInverse)
    XCTAssertEqual(product[0 * 4 + 0], 1, accuracy: 1e-5)
    XCTAssertEqual(product[0 * 4 + 1], 0, accuracy: 1e-5)
    XCTAssertEqual(product[0 * 4 + 2], 0, accuracy: 1e-5)
    XCTAssertEqual(product[0 * 4 + 3], 0, accuracy: 1e-5)
    
    XCTAssertEqual(product[1 * 4 + 1], 1, accuracy: 1e-5)
    XCTAssertEqual(product[1 * 4 + 2], 0, accuracy: 1e-5)
    XCTAssertEqual(product[1 * 4 + 3], 0, accuracy: 1e-5)
    XCTAssertEqual(product[2 * 4 + 2], 1, accuracy: 1e-5)
    XCTAssertEqual(product[2 * 4 + 3], 0, accuracy: 1e-5)
    XCTAssertEqual(product[3 * 4 + 3], 1, accuracy: 1e-5)
    
    // Invert the 4x4 matrix with Diagonalization.
    var diagonalizationDesc = DiagonalizationDescriptor()
    diagonalizationDesc.matrix = laplacian
    diagonalizationDesc.problemSize = 4
    let diagonalization = Diagonalization(descriptor: diagonalizationDesc)
    
    // Check that the eigenvalues match a previously evaluated result.
    XCTAssertEqual(diagonalization.eigenvalues[0], -3.618, accuracy: 1e-3)
    XCTAssertEqual(diagonalization.eigenvalues[1], -2.618, accuracy: 1e-3)
    XCTAssertEqual(diagonalization.eigenvalues[2], -1.382, accuracy: 1e-3)
    XCTAssertEqual(diagonalization.eigenvalues[3], -0.382, accuracy: 1e-3)
    
    // Materialize the eigenvalues into a square matrix, then invert it.
    let inverseΛ: [Float] = [
      1 / diagonalization.eigenvalues[0], 0, 0, 0,
      0, 1 / diagonalization.eigenvalues[1], 0, 0,
      0, 0, 1 / diagonalization.eigenvalues[2], 0,
      0, 0, 0, 1 / diagonalization.eigenvalues[3],
    ]
    
    // The eigenvectors are returned in column-major order.
    let Σ = transpose4x4(diagonalization.eigenvectors)
    var calculatedInverse = multiply4x4(inverseΛ, transpose4x4(Σ))
    calculatedInverse = multiply4x4(Σ, calculatedInverse)
    
    // Check the correctness of inversion by diagonalization.
    for elementID in 0..<16 {
      let expected = laplacianInverse[elementID]
      let actual = calculatedInverse[elementID]
      XCTAssertEqual(expected, actual, accuracy: 1e-5)
    }
  }
  
  // Directly invert the Laplacian for a 10-cell domain.
  func testLaplacianInverse10x10() throws {
    // The matrices are in row-major order.
    func multiply10x10(_ A: [Float], _ B: [Float]) -> [Float] {
      LinearAlgebraUtilities.matrixMultiply(matrixA: A, matrixB: B, n: 10)
    }
    func transpose10x10(_ A: [Float]) -> [Float] {
      var output = [Float](repeating: 0, count: 100)
      for m in 0..<10 {
        for n in 0..<10 {
          let value = A[n * 10 + m]
          output[m * 10 + n] = value
        }
      }
      return output
    }
    
    // Create a 10x10 matrix that represents an aperiodic Laplacian.
    var laplacian: [Float] = []
    laplacian += [-2, 1, 0, 0, 0, 0, 0, 0, 0, 0]
    laplacian += [1, -2, 1, 0, 0, 0, 0, 0, 0, 0]
    laplacian += [0, 1, -2, 1, 0, 0, 0, 0, 0, 0]
    laplacian += [0, 0, 1, -2, 1, 0, 0, 0, 0, 0]
    laplacian += [0, 0, 0, 1, -2, 1, 0, 0, 0, 0]
    laplacian += [0, 0, 0, 0, 1, -2, 1, 0, 0, 0]
    laplacian += [0, 0, 0, 0, 0, 1, -2, 1, 0, 0]
    laplacian += [0, 0, 0, 0, 0, 0, 1, -2, 1, 0]
    laplacian += [0, 0, 0, 0, 0, 0, 0, 1, -2, 1]
    laplacian += [0, 0, 0, 0, 0, 0, 0, 0, 1, -2]
    
    // Invert the 10x10 matrix with Diagonalization.
    var diagonalizationDesc = DiagonalizationDescriptor()
    diagonalizationDesc.matrix = laplacian
    diagonalizationDesc.problemSize = 10
    let diagonalization = Diagonalization(descriptor: diagonalizationDesc)
    
    // Check that the eigenvalues match a previously evaluated result.
    var expectedEigenvalues: [Float] = []
    expectedEigenvalues += [-3.919, -3.683, -3.310, -2.831, -2.285]
    expectedEigenvalues += [-1.715, -1.169, -0.690, -0.317, -0.081]
    for eigenpairID in 0..<10 {
      let expected = expectedEigenvalues[eigenpairID]
      let actual = diagonalization.eigenvalues[eigenpairID]
      XCTAssertEqual(expected, actual, accuracy: 1e-3)
    }
    
    // Materialize the eigenvalues into a square matrix, then invert it.
    var inverseΛ = [Float](repeating: .zero, count: 100)
    for diagonalID in 0..<10 {
      let eigenvalue = diagonalization.eigenvalues[diagonalID]
      let address = diagonalID * 10 + diagonalID
      inverseΛ[address] = 1 / eigenvalue
    }
    
    // The eigenvectors are returned in column-major order.
    let Σ = transpose10x10(diagonalization.eigenvectors)
    var laplacianInverse = multiply10x10(inverseΛ, transpose10x10(Σ))
    laplacianInverse = multiply10x10(Σ, laplacianInverse)
    
    // Evaluate the product of these two matrices.
    let product = multiply10x10(laplacian, laplacianInverse)
    for rowID in 0..<10 {
      let diagonal = product[rowID * 10 + rowID]
      XCTAssertEqual(diagonal, 1, accuracy: 1e-5)
      
      for columnID in (rowID + 1)..<10 {
        let entry = product[rowID * 10 + columnID]
        XCTAssertEqual(entry, 0, accuracy: 1e-5)
      }
    }
  }
  
  func testLaplacianPreconditioner() throws {
    let h: Float = 0.3
    
    // The matrices are in row-major order.
    func multiply10x10(_ A: [Float], _ B: [Float]) -> [Float] {
      LinearAlgebraUtilities.matrixMultiply(matrixA: A, matrixB: B, n: 10)
    }
    func transpose10x10(_ A: [Float]) -> [Float] {
      var output = [Float](repeating: 0, count: 100)
      for m in 0..<10 {
        for n in 0..<10 {
          let value = A[n * 10 + m]
          output[m * 10 + n] = value
        }
      }
      return output
    }
    
    // Create a 10x10 matrix that represents an aperiodic Laplacian.
    var laplacian: [Float] = []
    laplacian += [-2, 1, 0, 0, 0, 0, 0, 0, 0, 0]
    laplacian += [1, -2, 1, 0, 0, 0, 0, 0, 0, 0]
    laplacian += [0, 1, -2, 1, 0, 0, 0, 0, 0, 0]
    laplacian += [0, 0, 1, -2, 1, 0, 0, 0, 0, 0]
    laplacian += [0, 0, 0, 1, -2, 1, 0, 0, 0, 0]
    laplacian += [0, 0, 0, 0, 1, -2, 1, 0, 0, 0]
    laplacian += [0, 0, 0, 0, 0, 1, -2, 1, 0, 0]
    laplacian += [0, 0, 0, 0, 0, 0, 1, -2, 1, 0]
    laplacian += [0, 0, 0, 0, 0, 0, 0, 1, -2, 1]
    laplacian += [0, 0, 0, 0, 0, 0, 0, 0, 1, -2]
    
    // Adjust to match the 0.3 a.u. spacing.
    for entryID in laplacian.indices {
      var entry = laplacian[entryID]
      entry /= Float(h * h)
      laplacian[entryID] = entry
    }
    
    // Invert the 10x10 matrix with Diagonalization.
    var diagonalizationDesc = DiagonalizationDescriptor()
    diagonalizationDesc.matrix = laplacian
    diagonalizationDesc.problemSize = 10
    let diagonalization = Diagonalization(descriptor: diagonalizationDesc)
    
    // Materialize the eigenvalues into a square matrix, then invert it.
    var inverseΛ = [Float](repeating: .zero, count: 100)
    for diagonalID in 0..<10 {
      let eigenvalue = diagonalization.eigenvalues[diagonalID]
      let address = diagonalID * 10 + diagonalID
      inverseΛ[address] = 1 / eigenvalue
    }
    
    // The eigenvectors are returned in column-major order.
    let Σ = transpose10x10(diagonalization.eigenvectors)
    var laplacianInverse = multiply10x10(inverseΛ, transpose10x10(Σ))
    laplacianInverse = multiply10x10(Σ, laplacianInverse)
    
    // Visualize the exact inverse.
    XCTAssertEqual(laplacianInverse[0], -0.082, accuracy: 1e-3)
    XCTAssertEqual(laplacianInverse[1], -0.074, accuracy: 1e-3)
    XCTAssertEqual(laplacianInverse[2], -0.065, accuracy: 1e-3)
    
    // Visualize the expensive approximation to the inverse.
    var approximateInverse: [Float] = []
    for rowID in 0..<10 {
      let ri = (Float(rowID) + 0.5) * h
      for columnID in 0..<10 {
        let rj = (Float(columnID) + 0.5) * h
        let r = (ri - rj).magnitude
        let K = 1 / max(r, h)
        
        let diagonal = -2 / (h * h)
        let entry = (1 / diagonal) * K
        approximateInverse.append(entry)
      }
    }
    XCTAssertEqual(approximateInverse[0], -0.150, accuracy: 1e-3)
    XCTAssertEqual(approximateInverse[1], -0.150, accuracy: 1e-3)
    XCTAssertEqual(approximateInverse[2], -0.075, accuracy: 1e-3)
    
    // Visualize the efficient approximation to the inverse.
    var efficientInverse: [Float] = []
    for rowID in 0..<10 {
      let ri = (Float(rowID) + 0.5) * h
      for columnID in 0..<10 {
        let rj = (Float(columnID) + 0.5) * h
        let r = (ri - rj).magnitude
        
        // 0.8, 30, 0.2, 10 -> 0.8, 2.70, 0.2, 0.90
        // 0.6, 25, 0.4,  8 -> 0.6, 2.25, 0.4, 0.72
        var K: Float = .zero
        K += 0.6 * Float.exp(-2.25 * (r * r) / (h * h))
        K += 0.4 * Float.exp(-0.72 * (r * r) / (h * h))
        
        // Scale by (h * h) to visualize in comparison to the other matrices,
        // such as the exact inverse. Real-world usage will not modify the
        // diagonal, leaving the preconditioner close to the identity.
        let entry = (-h * h) * K
        efficientInverse.append(entry)
      }
    }
    XCTAssertEqual(efficientInverse[0], -0.090, accuracy: 1e-3)
    XCTAssertEqual(efficientInverse[1], -0.023, accuracy: 1e-3)
    XCTAssertEqual(efficientInverse[2], -0.002, accuracy: 1e-3)
  }
}
