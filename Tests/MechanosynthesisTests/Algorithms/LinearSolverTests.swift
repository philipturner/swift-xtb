import XCTest
import Mechanosynthesis
import Numerics

final class LinearSolverTests: XCTestCase {
  // Creates a 1D Laplacian operator with 2nd order accuracy.
  #if false
  // NOTE: Coded prematurely, just here as reference.
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
  #endif
  
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
  
  // Multigrid solver. There's currently a big unknown regarding how the grid
  // should treat domain boundaries.
  func testMultigrid() throws {
    
  }
}
