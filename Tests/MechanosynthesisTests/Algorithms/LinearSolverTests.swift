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
  
  // LAPACK solver as a reference implementation. This could be useful when
  // debugging the iterative solvers.
  func testDirectMatrixMethod() throws {
    // Define the matrix of equation coefficients.
    let coefficients: [Float] = [
      1, 1, 1,
      4, 2, 1,
      9, 3, 1,
    ]
    
    // Solve an equation where the first column is the RHS.
    do {
      let rightHandSide: [Float] = [1, 4, 9]
      let solution = LinearAlgebraUtilities
        .solveLinearSystem(matrix: coefficients, vector: rightHandSide, n: 3)
      XCTAssertEqual(solution[0], 1.000, accuracy: 1e-3)
      XCTAssertEqual(solution[1], 0.000, accuracy: 1e-3)
      XCTAssertEqual(solution[2], 0.000, accuracy: 1e-3)
    }
    
    // Solve an equation where the second column is the RHS.
    do {
      let rightHandSide: [Float] = [1, 2, 3]
      let solution = LinearAlgebraUtilities
        .solveLinearSystem(matrix: coefficients, vector: rightHandSide, n: 3)
      XCTAssertEqual(solution[0], 0.000, accuracy: 1e-3)
      XCTAssertEqual(solution[1], 1.000, accuracy: 1e-3)
      XCTAssertEqual(solution[2], 0.000, accuracy: 1e-3)
    }
    
    // Here, the RHS is a multiple of the third column.
    do {
      let rightHandSide: [Float] = [-2, -2, -2]
      let solution = LinearAlgebraUtilities
        .solveLinearSystem(matrix: coefficients, vector: rightHandSide, n: 3)
      XCTAssertEqual(solution[0], 0.000, accuracy: 1e-3)
      XCTAssertEqual(solution[1], 0.000, accuracy: 1e-3)
      XCTAssertEqual(solution[2], -2.000, accuracy: 1e-3)
    }
  }
  
  // Implementation of the algorithm from the INQ codebase, which chooses the
  // timestep based on the results of some integrals.
  func testSteepestDescent() throws {
    
  }
  
  // Implementation of weighted Jacobi, using a fixed timestep determined by
  // the grid spacing.
  func testWeightedJacobi() throws {
    
  }
  
  // Implementation of Gauss-Seidel, using a fixed timestep determined by the
  // grid spacing.
  //
  // This test does not cover the Gauss-Seidel red-black ordering scheme.
  // However, the results should reveal how one would go about coding GSRB.
  func testGaussSeidel() throws {
    
  }
  
  // Implementation of the algorithm from the INQ codebase, which chooses the
  // timestep based on the results of some integrals.
  func testConjugateGradient() throws {
    
  }
}
