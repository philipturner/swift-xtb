// Test the accuracy of Mehrstellen on various 3D functions.
//
// Hartree differential equation
//
// U (operator) = -v_{I}(r)
// ∇^2 (operator) Ψ(r) + 2(E - U (operator)) Ψ(r) = 0
// Given: v_{I}(r) = 1 / r
// Given: E = -1/2
// Solution: Ψ(r) = e^{-r} / (√π)
//
// Poisson differential equation
//
// ρ(r) = -Ψ(r) * Ψ(r)
// v_{H}(r) = ∫ ρ(r')dr' / |r - r'|
// ∇^2 (operator) v_{H}(r) = -4πρ(r)
// Given: Ψ(r) = e^{-r} / (√π)
// Solution: v_{H}(r) = e^{-2r}(1 + 1/r) - 1/r
//
// ======================================================================== //
// Methods
// ======================================================================== //
//
// Data for specific sample points:
// - Took the Laplacian in whichever differential equation is referenced.
// - Compared the Laplacian term to the other terms (the right hand side).
// - Subtracted RHS - LHS to create the residual.
// - Normalized the residual by dividing it by the RHS.
// - Reported the normalized residual to 7 decimal places.
//
// Data for integrals:
// - Encapsulated the hydrogen atom in a 16x16x16 Bohr uniform grid. The
//   finest grid spacing was 0.125 Bohr, and it took ~4 seconds to compute
//   all of the integrals.
// - Divided every Mehrstellen integral by ΨBΨ, to mirror the formula used
//   for Rayleigh quotients.
//
// All data is calculated in FP32, unless noted otherwise. Integrals are
// summed with FP32 integrands, and temporary FP64 accumulators.
//
// ======================================================================== //
// Results (Raw Data)
// ======================================================================== //
//
// Key:       2nd Order  4th Order  6th Order  8th Order  Mehrstellen
//
// Hartree potential, r ≈ 0.05
// h = 0.500, 0.3156697, 0.2372360, 0.2121769, 0.2002852, 0.0948272,
// h = 0.250, 0.1456027, 0.0889145, 0.0738506, 0.0672043, 0.0310335,
// h = 0.125, 0.0516573, 0.0203581, 0.0136607, 0.0109489, 0.0046150,
// h = 0.125, 0.0516985, 0.0203931, 0.0137204, 0.0110157, 0.0046477, (FP64)
//
// Hartree potential, r ≈ 0.33
// h = 0.500, 0.1402474, 0.0582366, 0.0352753, 0.0253011, 0.0101832,
// h = 0.250, 0.0386179, 0.0047414, 0.0006261, 0.0000823, 0.0053651,
// h = 0.125, 0.0105852, 0.0012434, 0.0009124, 0.0008343, 0.0008314,
// h = 0.125, 0.0105990, 0.0012587, 0.0009490, 0.0008495, 0.0008333, (FP64)
//
// Hartree potential, r ≈ 1
// h = 0.500, 0.0501792, 0.0297362, 0.0313445, 0.0327737, 0.0303154,
// h = 0.250, 0.0149480, 0.0031947, 0.0007120, 0.0001764, 0.0020250,
// h = 0.125, 0.0039315, 0.0003057, 0.0000095, 0.0000004, 0.0001522,
// h = 0.125, 0.0038992, 0.0002218, 0.0000145, 0.0000018, 0.0001264, (FP64)
//
// Hartree potential, r ≈ 2.5
// h = 0.500, 0.0596931, 0.0012162, 0.0005437, 0.0001839, 0.0005137,
// h = 0.250, 0.0149492, 0.0000688, 0.0000114, 0.0002007, 0.0000072,
// h = 0.125, 0.0034106, 0.0006545, 0.0000667, 0.0005778, 0.0001317,
// h = 0.125, 0.0037473, 0.0000040, 0.0000001, 0.0000000, 0.0000015, (FP64)
//
// Hartree potential, r ≈ 7
// h = 0.500, 4.8223810, 0.0777047, 0.7481260, 0.3463864, 0.1541908,
// h = 0.250, 4.2200661, 0.0791482, 3.2350295, 2.4064538, 0.2866600,
// h = 0.125, 10.243218, 8.0347290, 22.691066, 7.8276830, 2.8101072,
// h = 0.125, 0.2632410, 0.0004656, 0.0000007, 0.0000000, 0.0004544, (FP64)
//
// Wave function, r ≈ 0.05
// h = 0.500, 0.7630899, 0.7175785, 0.6992880, 0.6895115, 0.3422492,
// h = 0.250, 0.5242353, 0.4446174, 0.4138601, 0.3977126, 0.1917558,
// h = 0.125, 0.2038723, 0.0970867, 0.0590142, 0.0398512, 0.0038997,
// h = 0.125, 0.2038725, 0.0970849, 0.0590134, 0.0398527, 0.0038990, (FP64)
//
// Wave function, r ≈ 0.33
// h = 0.500, 0.1355058, 0.0080879, 0.0363433, 0.0581663, 0.0474944,
// h = 0.250, 0.0507169, 0.0224567, 0.0302432, 0.0390301, 0.0860758,
// h = 0.125, 0.0327208, 0.0267181, 0.0229551, 0.0177702, 0.0106808,
// h = 0.125, 0.0327156, 0.0267152, 0.0229603, 0.0177672, 0.0106801, (FP64)
//
// Wave function, r ≈ 1
// h = 0.500, 0.2025823, 0.2485496, 0.2730824, 0.2823939, 0.1626195,
// h = 0.250, 0.0600307, 0.0125117, 0.0031340, 0.0030094, 0.0071043,
// h = 0.125, 0.0156457, 0.0008615, 0.0000282, 0.0000040, 0.0004041,
// h = 0.125, 0.0156264, 0.0008276, 0.0000386, 0.0000249, 0.0003899, (FP64)
//
// Wave function, r ≈ 2.5
// h = 0.500, 0.0864438, 0.0010875, 0.0007733, 0.0007574, 0.0007132,
// h = 0.250, 0.0216032, 0.0000785, 0.0000801, 0.0000197, 0.0000681,
// h = 0.125, 0.0054886, 0.0000284, 0.0000383, 0.0003043, 0.0000942,
// h = 0.125, 0.0053842, 0.0000046, 0.0000000, 0.0000000, 0.0000040, (FP64)
//
// Wave function, r ≈ 7
// h = 0.500, 0.0062805, 0.0000452, 0.0000197, 0.0000140, 0.0001070,
// h = 0.250, 0.0015472, 0.0000248, 0.0000551, 0.0000555, 0.0000063,
// h = 0.125, 0.0003074, 0.0001595, 0.0001967, 0.0001691, 0.0000276,
// h = 0.125, 0.0003937, 0.0000002, 0.0000000, 0.0000000, 0.0000004, (FP64)
//
// Wave function, r ≈ 12
// h = 0.500, 0.0094234, 0.0001468, 0.0000252, 0.0000232, 0.0000074,
// h = 0.250, 0.0022993, 0.0000669, 0.0000624, 0.0000746, 0.0000263,
// h = 0.125, 0.0002278, 0.0003585, 0.0005138, 0.0005083, 0.0001787,
// h = 0.125, 0.0005881, 0.0000005, 0.0000000, 0.0000000, 0.0000001, (FP64)
//
// Euclidean Norm
// h = 0.500, 0.9782071, 0.9782071, 0.9782071, 0.9782071, 1.0000000,
// h = 0.250, 0.9947214, 0.9947214, 0.9947214, 0.9947214, 1.0000000,
// h = 0.125, 0.9986905, 0.9986905, 0.9986905, 0.9986905, 1.0000000,
// h = 0.125, 0.9986905, 0.9986905, 0.9986905, 0.9986905, 1.0000000, (FP64)
//
// Kinetic Energy
// h = 0.500, 0.4670562, 0.4825225, 0.4845980, 0.4851608, 0.4658559,
// h = 0.250, 0.4918282, 0.4973585, 0.4977855, 0.4978775, 0.4905346,
// h = 0.125, 0.4979922, 0.4996298, 0.4997187, 0.4997239, 0.4975142,
// h = 0.125, 0.4979936, 0.4996396, 0.4997066, 0.4997194, 0.4975183, (FP64)
//
// Potential Energy
// h = 0.500, -0.9385867, -0.9385867, -0.9385867, -0.9385867, -0.9239331,
// h = 0.250, -0.9841995, -0.9841995, -0.9841995, -0.9841995, -0.9776164,
// h = 0.125, -0.9960209, -0.9960209, -0.9960209, -0.9960209, -0.9939297,
// h = 0.125, -0.9960208, -0.9960208, -0.9960208, -0.9960208, -0.9939297, (FP64)
//
// Rayleigh Quotient
// h = 0.500, -0.4715305, -0.4560642, -0.4539886, -0.4534259, -0.4580772,
// h = 0.250, -0.4923713, -0.4868410, -0.4864140, -0.4863219, -0.4870818,
// h = 0.125, -0.4980287, -0.4963911, -0.4963022, -0.4962970, -0.4964156,
// h = 0.125, -0.4980272, -0.4963813, -0.4963142, -0.4963014, -0.4964114, (FP64)
//
// Atomic Radius
// h = 0.500, 1.0006634, 1.0006634, 1.0006634, 1.0006634, 1.0498127,
// h = 0.250, 1.0000311, 1.0000311, 1.0000311, 1.0000311, 1.0122266,
// h = 0.125, 0.9999850, 0.9999850, 0.9999850, 0.9999850, 1.0030289,
// h = 0.125, 0.9999849, 0.9999849, 0.9999849, 0.9999849, 1.0030289, (FP64)
func testMehrstellenDiscretization() throws {
  typealias Real = Float
  
  // The analytical solutions to the differential equations.
  func waveFunction(r: SIMD3<Real>) -> Real {
    let radius = (r * r).sum().squareRoot()
    return Real.exp(-radius) / Real.pi.squareRoot()
  }
  func ionicPotential(r: SIMD3<Real>) -> Real {
    let radius = (r * r).sum().squareRoot()
    return 1 / radius
  }
  func chargeDensity(r: SIMD3<Real>) -> Real {
    let waveFunctionValue = waveFunction(r: r)
    return -waveFunctionValue * waveFunctionValue
  }
  func hartreePotential(r: SIMD3<Real>) -> Real {
    let radius = (r * r).sum().squareRoot()
    var output: Real = .zero
    output += Real.exp(-2 * radius) * (1 + 1 / radius)
    output += -1 / radius
    return output
  }
  
  // The available techniques for sampling discretized functions.
  enum FiniteDifferenceMethod {
    case secondOrderLaplacian
    case fourthOrderLaplacian
    case sixthOrderLaplacian
    case eighthOrderLaplacian
    case mehrstellenLeftHandSide
    case mehrstellenRightHandSide
  }
  
  // A self-contained function for executing finite differences.
  func sample(
    method finiteDifferenceMethod: FiniteDifferenceMethod,
    spacing: Real,
    position: SIMD3<Real>,
    function: (SIMD3<Real>) -> Real
  ) -> Real {
    var accumulator: Real = .zero
    
    // Accumulate the central point.
    let centerValue = function(position)
    switch finiteDifferenceMethod {
    case .secondOrderLaplacian:
      accumulator += 3 * Real(-2) * centerValue
    case .fourthOrderLaplacian:
      accumulator += 3 * Real(-5.0 / 2) * centerValue
    case .sixthOrderLaplacian:
      accumulator += 3 * Real(-49.0 / 18) * centerValue
    case .eighthOrderLaplacian:
      accumulator += 3 * Real(-205.0 / 72) * centerValue
    case .mehrstellenLeftHandSide:
      accumulator += Real(-4) * centerValue
    case .mehrstellenRightHandSide:
      accumulator += Real(1.0 / 2) * centerValue
    }
    
    // Iterate over the faces.
    for faceID in 0..<6 {
      let coordinateID = faceID / 2
      var coordinateShift: SIMD3<Real> = .zero
      coordinateShift[coordinateID] = (faceID % 2 == 0) ? -1 : 1
      
      // Accumulate the central point, plus h.
      let firstNeighborPosition = position + coordinateShift * spacing
      let firstNeighborValue = function(firstNeighborPosition)
      switch finiteDifferenceMethod {
      case .secondOrderLaplacian:
        accumulator += Real(1) * firstNeighborValue
      case .fourthOrderLaplacian:
        accumulator += Real(4.0 / 3) * firstNeighborValue
      case .sixthOrderLaplacian:
        accumulator += Real(3.0 / 2) * firstNeighborValue
      case .eighthOrderLaplacian:
        accumulator += Real(8.0 / 5) * firstNeighborValue
      case .mehrstellenLeftHandSide:
        accumulator += Real(1.0 / 3) * firstNeighborValue
      case .mehrstellenRightHandSide:
        accumulator += Real(1.0 / 12) * firstNeighborValue
      }
      
      // Accumulate the central point, plus 2 * h.
      let secondNeighborPosition = position + 2 * coordinateShift * spacing
      let secondNeighborValue = function(secondNeighborPosition)
      switch finiteDifferenceMethod {
      case .fourthOrderLaplacian:
        accumulator += Real(-1.0 / 12) * secondNeighborValue
      case .sixthOrderLaplacian:
        accumulator += Real(-3.0 / 20) * secondNeighborValue
      case .eighthOrderLaplacian:
        accumulator += Real(-1.0 / 5) * secondNeighborValue
      default:
        break
      }
      
      // Accumulate the central point, plus 3 * h.
      let thirdNeighborPosition = position + 3 * coordinateShift * spacing
      let thirdNeighborValue = function(thirdNeighborPosition)
      switch finiteDifferenceMethod {
      case .sixthOrderLaplacian:
        accumulator += Real(1.0 / 90) * thirdNeighborValue
      case .eighthOrderLaplacian:
        accumulator += Real(8.0 / 315) * thirdNeighborValue
      default:
        break
      }
      
      // Accumulate the central point, plus 4 * h.
      let fourthNeighborPosition = position + 4 * coordinateShift * spacing
      let fourthNeighborValue = function(fourthNeighborPosition)
      switch finiteDifferenceMethod {
      case .eighthOrderLaplacian:
        accumulator += Real(-1.0 / 560) * fourthNeighborValue
      default:
        break
      }
    }
    
    // Iterate over the edges.
    if finiteDifferenceMethod == .mehrstellenLeftHandSide {
      func createEdgePermutations() -> [SIMD3<Real>] {
        var permutations: [SIMD3<Real>] = []
        permutations.append(SIMD3(1, 1, 0))
        permutations.append(SIMD3(1, 0, 1))
        permutations.append(SIMD3(0, 1, 1))
        permutations.append(SIMD3(1, -1, 0))
        permutations.append(SIMD3(1, 0, -1))
        permutations.append(SIMD3(0, 1, -1))
        permutations += permutations.map(-)
        return permutations
      }
      let edgePermutations = createEdgePermutations()
      
      for coordinateShift in edgePermutations {
        let edgePosition = position + coordinateShift * spacing
        let edgeValue = function(edgePosition)
        accumulator += Real(1.0 / 6) * edgeValue
      }
    }
    
    // Scale by the grid spacing.
    if finiteDifferenceMethod != .mehrstellenRightHandSide {
      accumulator /= Real(spacing * spacing)
    }
    
    return accumulator
  }
  
  // Specific locations within the domain to report the residual at.
  var samplePoints: [SIMD3<Real>] = []
  samplePoints.append(SIMD3(0.046279, 0.017622, -0.019439))// r ≈ 0.05
  samplePoints.append(SIMD3(0.034873, -0.082556, -0.249546)) // r ≈ 0.33
  samplePoints.append(SIMD3(-0.064516, 0.126008, -0.800164)) // r ≈ 1
  samplePoints.append(SIMD3(-2.054362, 0.571690, 1.3407375)) // r ≈ 2.5
  samplePoints.append(SIMD3(-1.287692, 6.581067, -4.103087)) // r ≈ 7
  samplePoints.append(SIMD3(2.442969, -5.903516, -10.410264)) // r ≈ 12
  
  // Query various expectation values of the hydrogen electron.
  let resolutions: [Real] = [
    // Deactivate the unit test entirely, as the latency is not acceptable.
    // 1.0 / 2,
    // 1.0 / 4,
    // 1.0 / 8,
  ]
  
  // Iterate over the computationally feasible resolutions.
  for h in resolutions {
    // Iterate over the difference methods.
    var expectationValues: [Real] = []
    let finiteDifferenceMethods: [FiniteDifferenceMethod] = [
      .secondOrderLaplacian,
      .fourthOrderLaplacian,
      .sixthOrderLaplacian,
      .eighthOrderLaplacian,
      .mehrstellenLeftHandSide
    ]
    for method in finiteDifferenceMethods {
      // Encapsulate the hydrogen atom in a 16x16x16 Bohr simulation box.
      let gridSize = Int32(Float(16 / h).rounded(.toNearestOrEven))
      
      // Iterate over the cells.
      var accumulator: Double = .zero
      var mehrstellenCorrection: Double = .zero
      for indexZ in 0..<gridSize {
        for indexY in 0..<gridSize {
          for indexX in 0..<gridSize {
            // Locate the cell.
            let cellIndices = SIMD3(indexX, indexY, indexZ)
            var cellPosition = SIMD3<Real>(cellIndices)
            cellPosition = h * (cellPosition + 0.5)
            
            // Generate the position, relative to the nucleus.
            var nucleusPosition = SIMD3(repeating: Real(gridSize))
            nucleusPosition = h * (nucleusPosition * 0.5)
            let r = SIMD3<Real>(cellPosition - nucleusPosition)
            
            // Compute the integrand at this point in real space.
            let Ψ = waveFunction(r: r)
            var rΨ: Real
            if method == .mehrstellenLeftHandSide {
              rΨ = sample(
                method: .mehrstellenRightHandSide,
                spacing: h,
                position: r
              ) { r in
                // Be careful to not reference the Ψ from the outer scope.
                let radius = (r * r).sum().squareRoot()
                return Real(2.0 / 3) * radius * Ψ
              }
            } else {
              let radius = (r * r).sum().squareRoot()
              rΨ = Real(2.0 / 3) * radius * Ψ
            }
            let integrand = Ψ * rΨ
            
            // Accumulate the integral.
            let BΨ = sample(
              method: .mehrstellenRightHandSide,
              spacing: h,
              position: r,
              function: waveFunction(r:))
            let drTerm = Real(h * h * h)
            accumulator += Double(integrand * drTerm)
            mehrstellenCorrection += Double(Ψ * BΨ * drTerm)
          }
        }
      }
      
      // Add the integral to the list.
      if method == .mehrstellenLeftHandSide {
        let expectationValue = Real(accumulator / mehrstellenCorrection)
        expectationValues.append(expectationValue)
      } else {
        let expectationValue = Real(accumulator)
        expectationValues.append(expectationValue)
      }
    }
  }
}
