import XCTest
import Mechanosynthesis
import Numerics

// Test an alternative to Cholesky decomposition that's more parallel.
final class OrthogonalizeTests: XCTestCase {
  typealias Real = Float
  
  // Use real numbers instead of complex numbers.
  func testDimension3() throws {
    // TODO: Formulate this into a function that handles 3, 100, and 37. We
    // may need to reduce basis size and/or electron count to achieve reasonable
    // test latency.
    let numPoints = 3
    let spacing: Real = 0.44428829
    let gridWidth = Int((6.3 / spacing).rounded(.up))
    
    var phi: [[Real]] = []
    for _ in 0..<numPoints {
      var phiVector: [Real] = []
      for _ in 0..<(gridWidth * gridWidth * gridWidth) {
        let value = Real.random(in: 0..<1)
        phiVector.append(value)
      }
      
      var norm: Real = .zero
      for cellID in 0..<(gridWidth * gridWidth * gridWidth) {
        let d3r = spacing * spacing * spacing
        let Ψ = phiVector[cellID]
        norm += Ψ * 1 * Ψ * d3r
      }
      
      let normalizationFactor = 1 / norm.squareRoot()
      for cellID in 0..<(gridWidth * gridWidth * gridWidth) {
        phiVector[cellID] *= normalizationFactor
      }
      phi.append(phiVector)
    }
  }
}
