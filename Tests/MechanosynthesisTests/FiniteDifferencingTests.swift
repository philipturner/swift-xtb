import XCTest
import Mechanosynthesis
import Numerics

// The purpose of this test file is to gain experience with asymmetric finite
// differencing and higher-order finite differencing. This is a key component
// of the variable-resolution basis method. Examine the tradeoff between FD
// quality, grid resolution, and code complexity.
final class FiniteDifferencingTests: XCTestCase {
  // Test various orders of finite differencing with a 1D function on a uniform
  // grid. Demonstrate an improvement in quality, like the multipole expansion.
  //
  // Graph of the test function and its derivatives:
  // https://www.desmos.com/calculator/davk2k5ahe
  func testSymmetricDifference() throws {
    func Ψ(x: Float) -> Float {
      let expTerm = Float.exp(-x * x)
      let trigTerm = Float.sin(x)
      return expTerm * trigTerm
    }
    
    func dΨ_dx(x: Float) -> Float {
      let expTerm = Float.exp(-x * x)
      let trigTerm = -2 * x * Float.sin(x) + Float.cos(x)
      return expTerm * trigTerm
    }
    
    func d2Ψ_dx2(x: Float) -> Float {
      let expTerm = Float.exp(-x * x)
      let trigTerm = (4 * x * x - 3) * Float.sin(x) - 4 * x * Float.cos(x)
      return expTerm * trigTerm
    }
    
    var samplePoints: [Float] = []
    samplePoints.append(0)
    samplePoints.append(0.1)
    samplePoints.append(1.5)
    samplePoints.append(-2)
    
    // These are reference points to benchmark each finite difference against.
    var expectedSamples: [SIMD3<Float>] = []
    expectedSamples.append(SIMD3(0.0, 1.0, -0.0))
    expectedSamples.append(SIMD3(0.098840065, 0.9653357, -0.6866081))
    expectedSamples.append(SIMD3(0.105135195, -0.30794996, 0.58607733))
    expectedSamples.append(SIMD3(-0.016654363, -0.074239455, -0.2774827))
    
    for pointID in samplePoints.indices {
      let point = samplePoints[pointID]
      let zerothDerivative = Ψ(x: point)
      let firstDerivative = dΨ_dx(x: point)
      let secondDerivative = d2Ψ_dx2(x: point)
      
      let expected = expectedSamples[pointID]
      XCTAssertEqual(zerothDerivative, expected[0], accuracy: 1e-3)
      XCTAssertEqual(firstDerivative, expected[1], accuracy: 1e-3)
      XCTAssertEqual(secondDerivative, expected[2], accuracy: 1e-3)
    }
    
    
  }
  
  // Test asymmetric 2nd-order FD with a 1D variable-resolution grid. Test
  // some edge cases that jump multiple resolution levels.
  
  // Evaluate the kinetic energy of an N atom ansatz, with uniform and
  // variable-resolution grids. Use only two resolution levels to make
  // the code simpler.
}
