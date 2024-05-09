import XCTest
import Mechanosynthesis
import Numerics

final class FiniteDifferencingTests: XCTestCase {
  // Test various orders of finite differencing with a 1D function on a uniform
  // grid. Demonstrate an improvement in quality, like the multipole expansion.
  //
  // Graph of the test function and its derivatives:
  // https://www.desmos.com/calculator/davk2k5ahe
  func testSymmetricDifference() throws {
    typealias Real = Float
    
    // Program settings.
    let h: Real = 0.05 // must divide evenly into 0.1
    
    // It looks like 4th-order FD with a large grid spacing is optimal for FP32.
    // Mehrstellen discretization might be interesting to explore as well.
    
    /*
     FP32
     
     h = 0.1
     [0.0, 7.45058e-07, 9.934107e-07, 1.0844733e-06, 1.1164044e-06, 0.0]
     [0.6866081, 0.00653404, 0.0001347065, 1.4305115e-06, 2.9802322e-06, 0.0]
     [0.58607733, 0.00504905, 2.2351742e-05, 1.9669533e-06, 5.364418e-07, 0.0]
     [0.2774827, 0.0008057058, 3.72231e-05, 2.0563602e-06, 2.7418137e-06, 0.0]
     
     h = 0.05
     [0.0, 2.980232e-06, 3.72529e-06, 4.0233135e-06, 4.1723247e-06, 0.0]
     [0.6866081, 0.0016419888, 1.1205673e-05, 3.3974648e-06, 3.874302e-06, 0.0]
     [0.58607733, 0.001242578, 2.6285648e-05, 3.0994415e-05, 3.2901764e-05, 0.0]
     [0.2774827, 0.00020304322, 1.2218952e-06, 1.66893e-06, 1.937151e-06, 0.0]
     
     h = 0.025
     [0.0, 1.1920928e-05, 1.490116e-05, 1.6358163e-05, 1.7294806e-05, 0.0]
     [0.6866081, 0.0003798604, 4.0352345e-05, 4.6491623e-05, 5.0008297e-05, 0.0]
     [0.58607733, 0.0003067851, 1.3113022e-05, 1.5735626e-05, 1.758337e-05, 0.0]
     [0.2774827, 6.631017e-05, 2.2232533e-05, 2.5868416e-05, 2.8282404e-05, 0.0]
     
     h = 0.0125
     [0.0, 4.7683712e-05, 6.3578285e-05, 7.23203e-05, 7.815304e-05, 0.0]
     [0.6866081, 1.0311604e-05, 0.0001128912, 0.00012135506, 0.00012606382, 0.0]
     [0.58607733, 0.00023525953, 0.00021141768, 0.00024056435, 0.00025856495, 0.0]
     [0.2774827, 2.4586916e-05, 2.4586916e-05, 3.2544136e-05, 3.7670135e-05, 0.0]
     
     h = 0.00625
     [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
     [0.6866081, 0.00022810698, 0.00030761957, 0.00033086538, 0.00034421682, 0.0]
     [0.58607733, 0.0008137822, 0.0010045767, 0.0010660291, 0.0010953546, 0.0]
     [0.2774827, 0.0004403293, 0.0005515516, 0.00059849024, 0.00062435865, 0.0]
     */
    
    /*
     FP64
     
     h = 0.1
     [0.0, 0.0, 0.0, 6.167906e-17, 1.409807e-16, 0.0]
     [0.6866081, 0.006536034, 0.00013736018, 4.413789e-06, 1.893477e-07, 0.0]
     [0.58607733, 0.00504985, 2.2994873e-05, 1.5228054e-06, 1.1035841e-07, 0.0]
     [0.2774827, 0.000807538, 3.957382e-05, 5.298252e-07, 1.9880988e-08, 0.0]
     
     h = 0.05
     [0.0, 5.551115e-15, 7.4014865e-15, 8.450031e-15, 9.163745e-15, 0.0]
     [0.6866081, 0.0016405355, 8.702746e-06, 7.098031e-08, 7.7709505e-10, 0.0]
     [0.58607733, 0.0012635094, 1.395899e-06, 2.4966019e-08, 4.5177173e-10, 0.0]
     [0.2774827, 0.00020374965, 2.4868575e-06, 8.058645e-09, 8.583023e-11, 0.0]
     
     h = 0.025
     [0.0, 2.220446e-14, 2.9605946e-14, 3.3553408e-14, 3.605141e-14, 0.0]
     [0.6866081, 0.00041054323, 5.4578123e-07, 1.1171124e-09, 3.075984e-12, 0.0]
     [0.58607733, 0.00031594228, 8.658697e-08, 3.9484593e-10, 1.8708368e-12, 0.0]
     [0.2774827, 5.1054136e-05, 1.5563728e-07, 1.2509732e-10, 2.624012e-13, 0.0]
     
     h = 0.0125
     [0.0, 8.881784e-14, 1.110223e-13, 1.2039752e-13, 1.2523492e-13, 0.0]
     [0.6866081, 0.00010266141, 3.4140303e-08, 1.7305823e-11, 1.8185453e-13, 0.0]
     [0.58607733, 7.898962e-05, 5.4018696e-09, 5.655809e-12, 5.4933835e-13, 0.0]
     [0.2774827, 1.2770832e-05, 9.730356e-09, 1.7041923e-12, 2.6167957e-13, 0.0]
     
     h = 0.00625
     [0.0, 3.5527137e-13, 4.7369514e-13, 5.3290705e-13, 5.6779977e-13, 0.0]
     [0.6866081, 2.5666954e-05, 2.1337345e-09, 2.6623148e-13, 5.6121774e-13, 0.0]
     [0.58607733, 1.9747657e-05, 3.3547043e-10, 2.349232e-12, 2.4406033e-12, 0.0]
     [0.2774827, 3.1931638e-06, 6.0775907e-10, 5.0498494e-13, 5.8558713e-13, 0.0]
     */
    
    func Ψ(x: Real) -> Real {
      let expTerm = Real.exp(-x * x)
      let trigTerm = Real.sin(x)
      return expTerm * trigTerm
    }
    
    func dΨ_dx(x: Real) -> Real {
      let expTerm = Real.exp(-x * x)
      let trigTerm = -2 * x * Real.sin(x) + Real.cos(x)
      return expTerm * trigTerm
    }
    
    func d2Ψ_dx2(x: Real) -> Real {
      let expTerm = Real.exp(-x * x)
      let trigTerm = (4 * x * x - 3) * Real.sin(x) - 4 * x * Real.cos(x)
      return expTerm * trigTerm
    }
    
    var samplePoints: [Real] = []
    samplePoints.append(0)
    samplePoints.append(0.1)
    samplePoints.append(1.5)
    samplePoints.append(-2)
    
    // These are reference points to benchmark each finite difference against.
    var expectedSamples: [SIMD3<Real>] = []
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
    
    // Set the origin so that each sample point coincides with the center of a
    // cell.
    let origin: Real = 0.1 - h / 2
    
    // Test the various finite differences.
    for pointID in samplePoints.indices {
      let pointX = samplePoints[pointID]
      let numCells = (pointX - origin) / h
      
      func createFiniteDifference(cellOffset: Int) -> Real {
        let neighborCellLeft = numCells - Real(cellOffset)
        let neighborCellRight = numCells + Real(cellOffset)
        let neighborXLeft = origin + neighborCellLeft * h
        let neighborXMiddle = origin + numCells * h
        let neighborXRight = origin + neighborCellRight * h
        
        var dΨdxLeft = Ψ(x: neighborXMiddle) - Ψ(x: neighborXLeft)
        var dΨdxRight = Ψ(x: neighborXRight) - Ψ(x: neighborXMiddle)
        if Float.random(in: 0..<1) > 2 {
          // Suppress compiler optimizations.
          dΨdxLeft = .nan
          dΨdxRight = .nan
        }
        return dΨdxRight - dΨdxLeft
      }
      
      var estimates: [Real] = []
      estimates.append(.zero)
      
      let difference2 = createFiniteDifference(cellOffset: 1)
      let difference4 = createFiniteDifference(cellOffset: 2)
      let difference6 = createFiniteDifference(cellOffset: 3)
      let difference8 = createFiniteDifference(cellOffset: 4)
      
      typealias Accumulator = Double
      
      // Symmetric second-order.
      do {
        var result: Accumulator = .zero
        result += Accumulator(1 * difference2)
        estimates.append(Real(result) / (h * h))
      }
      
      // Symmetric fourth-order.
      do {
        var result: Accumulator = .zero
        result += Accumulator(4 / 3 * difference2)
        result += Accumulator(-1 / 12 * difference4)
        estimates.append(Real(result) / (h * h))
      }
      
      // Symmetric sixth-order.
      do {
        var result: Accumulator = .zero
        result += Accumulator(3 / 2 * difference2)
        result += Accumulator(-3 / 20 * difference4)
        result += Accumulator(1 / 90 * difference6)
        estimates.append(Real(result) / (h * h))
      }
      
      // Symmetric eighth-order.
      do {
        var result: Accumulator = .zero
        result += Accumulator(8 / 5 * difference2)
        result += Accumulator(-1 / 5 * difference4)
        result += Accumulator(8 / 315 * difference6)
        result += Accumulator(-1 / 560 * difference8)
        estimates.append(Real(result) / (h * h))
      }
      
      estimates.append(d2Ψ_dx2(x: pointX))
      
      let errors = estimates.map { Float($0 - d2Ψ_dx2(x: pointX)).magnitude }
      
      XCTAssertLessThan(errors[0], 1)
      XCTAssertLessThan(errors[1], 0.002)
      XCTAssertLessThan(errors[2], 3e-5)
      XCTAssertLessThan(errors[3], 4e-5)
      XCTAssertLessThan(errors[4], 4e-5)
      XCTAssertEqual(errors[5], 0)
    }
  }
  
  // Test formulae for quadratic interpolation (2nd order FVM).
  //
  // f(x) = ax^2 + bx + c
  //
  // | x0^2 x0 1 | | a | = | y0 |
  // | x1^2 x1 1 | | b | = | y1 |
  // | x2^2 x2 1 | | c | = | y2 |
  //
  // (x0, y0) = ( -1, y0)
  // (x1, y1) = (  0, y1)
  // (x2, y2) = (3/2, y2)
  //
  // | a | = |  2/5 -2/3 4/15 | | y0 |
  // | b | = | -3/5  1/3 4/15 | | y1 |
  // | c | = |    0    1    0 | | y2 |
  //
  // (x0, y0) = (-1, y0)
  // (x1, y1) = ( 0, y1)
  // (x2, y2) = ( 1, y2)
  //
  // | a | = |  1/2 -1 1/2 | | y0 |
  // | b | = | -1/2  0 1/2 | | y1 |
  // | c | = |    0  1   0 | | y2 |
  func testQuadraticInterpolation() throws {
    
  }
}
