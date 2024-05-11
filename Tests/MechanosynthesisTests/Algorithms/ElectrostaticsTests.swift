import XCTest
import Mechanosynthesis
import Numerics

final class ElectrostaticsTests: XCTestCase {
  // This test defines the origin to the "center of charge", to have defined
  // behavior when the quadrupole depends on the origin:
  // https://phys.libretexts.org/Bookshelves/Mathematical_Physics_and_Pedagogy/Mathematical_Methods/The_Multipole_Expansion
  //
  // Also reproducing the octupole expansion. This has the best tradeoff
  // between computational overhead and improvement in quality:
  // https://physics.stackexchange.com/questions/269753/electric-octupole-moment-in-cartesian-coordinates
  //
  // For reference, the hexadecapole expansion is on this page (formula 19):
  //
  // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7017380/
  func testMultipoleExpansion() throws {
    // Create a list of point charges, compute the multipole expansion, and
    // compare to results from direct integration. Show how the accuracy
    // improves from cutoff -> monopole -> dipole -> quadrupole -> octupole.
    
    // Program settings.
    let ionCharge: Int = 1
    let samplePosition: SIMD3<Float> = [-1.5, 0, 0]
    
    // Create the list of charges.
    var pointCharges: [SIMD4<Float>] = []
    pointCharges.append(SIMD4(-0.1, 0.2, 0.2, 7))
    for _ in 0..<2 {
      pointCharges.append(SIMD4(0.1, 0.2, 0, -1))
      pointCharges.append(SIMD4(-0.3, 0, -0.5, -1))
      pointCharges.append(SIMD4(-0.01, 0.4, 0.03, -1))
    }
    for _ in 0..<(1 - ionCharge) {
      pointCharges.append(SIMD4(0.3, 0.4, 0.5, -1))
    }
    
    // Find the origin of the multipole expansion.
    var origin: SIMD3<Float> = .zero
    var totalMass: Float = .zero
    for charge in pointCharges {
      let position = unsafeBitCast(charge, to: SIMD3<Float>.self)
      let mass = charge.w.magnitude
      origin += position * mass
      totalMass += mass
    }
    origin /= totalMass
    
    // Pre-compute some quantities.
    let r = samplePosition - origin
    let rNorm = (r * r).sum().squareRoot()
    let rHat = r / rNorm
    
    // Record each estimate to the electric potential, starting with 'cutoff'.
    var estimates: [Float] = []
    var potential: Float = .zero
    estimates.append(potential)
    
    // Compute the monopole expansion.
    do {
      var totalCharge: Float = .zero
      for charge in pointCharges {
        totalCharge += charge.w
      }
      potential += totalCharge / rNorm
    }
    estimates.append(potential)
    
    // Compute the dipole expansion.
    do {
      var dipoleMoment: SIMD3<Float> = .zero
      for charge in pointCharges {
        let position = unsafeBitCast(charge, to: SIMD3<Float>.self)
        let rPrime = position - origin
        dipoleMoment += charge.w * rPrime
      }
      let dotProduct = (rHat * dipoleMoment).sum()
      potential += dotProduct / (rNorm * rNorm)
    }
    estimates.append(potential)
    
    // Compute the quadrupole expansion.
    do {
      var quadrupoleMoment: (SIMD3<Float>, SIMD3<Float>, SIMD3<Float>)
      quadrupoleMoment = (.zero, .zero, .zero)
      
      for charge in pointCharges {
        let position = unsafeBitCast(charge, to: SIMD3<Float>.self)
        let rPrime = position - origin
        let r2 = (rPrime * rPrime).sum()
        
        var matrix: (SIMD3<Float>, SIMD3<Float>, SIMD3<Float>)
        matrix = (
          3 * rPrime.x * rPrime,
          3 * rPrime.y * rPrime,
          3 * rPrime.z * rPrime)
        matrix.0[0] -= r2
        matrix.1[1] -= r2
        matrix.2[2] -= r2
        
        quadrupoleMoment.0 += charge.w * matrix.0
        quadrupoleMoment.1 += charge.w * matrix.1
        quadrupoleMoment.2 += charge.w * matrix.2
      }
      
      let M_rHat = SIMD3(
        (quadrupoleMoment.0 * rHat).sum(),
        (quadrupoleMoment.1 * rHat).sum(),
        (quadrupoleMoment.2 * rHat).sum())
      let matrixProduct = 1 / 2 * (rHat * M_rHat).sum()
      potential += matrixProduct / (rNorm * rNorm * rNorm)
    }
    estimates.append(potential)
    
    // Compute the octupole expansion.
    do {
      var octupoleMoment: (
        (SIMD3<Float>, SIMD3<Float>, Float),
        (Void, SIMD3<Float>, Float),
        (Void, Void, Float))
      octupoleMoment = (
        (.zero, .zero, .zero),
        ((), .zero, .zero),
        ((), (), .zero))
      
      for charge in pointCharges {
        let position = unsafeBitCast(charge, to: SIMD3<Float>.self)
        let rPrime = position - origin
        let r2 = (rPrime * rPrime).sum()
        
        var tensor: (
          (SIMD3<Float>, SIMD3<Float>, Float),
          (Void, SIMD3<Float>, Float),
          (Void, Void, Float))
        tensor = (
          (15 * rPrime.x * rPrime.x * rPrime,
           15 * rPrime.x * rPrime.y * rPrime,
           15 * rPrime.x * rPrime.z * rPrime.z),
          ((),
           15 * rPrime.y * rPrime.y * rPrime,
           15 * rPrime.y * rPrime.z * rPrime.z),
          ((),
           (),
           15 * rPrime.z * rPrime.z * rPrime.z)
        )
        
        tensor.0.0 -= 3 * rPrime.x * r2
        tensor.0.0[0] -= 3 * rPrime.x * r2
        tensor.0.1[0] -= 3 * rPrime.x * r2
        // tensor.0.2[0] -= 3 * rPrime.x * r2
        tensor.0.0[0] -= 3 * rPrime.x * r2
        // tensor.1.0[0] -= 3 * rPrime.x * r2
        // tensor.2.0[0] -= 3 * rPrime.x * r2
        
        tensor.1.1 -= 3 * rPrime.y * r2
        // tensor.1.0[1] -= 3 * rPrime.y * r2
        tensor.1.1[1] -= 3 * rPrime.y * r2
        // tensor.1.2[1] -= 3 * rPrime.y * r2
        tensor.0.1[1] -= 3 * rPrime.y * r2
        tensor.1.1[1] -= 3 * rPrime.y * r2
        // tensor.2.1[1] -= 3 * rPrime.y * r2
        
        tensor.2.2 -= 3 * rPrime.z * r2
        // tensor.2.0[2] -= 3 * rPrime.z * r2
        // tensor.2.1[2] -= 3 * rPrime.z * r2
        tensor.2.2 -= 3 * rPrime.z * r2
        tensor.0.2 -= 3 * rPrime.z * r2
        tensor.1.2 -= 3 * rPrime.z * r2
        tensor.2.2 -= 3 * rPrime.z * r2
        
        octupoleMoment.0.0 += charge.w * tensor.0.0
        octupoleMoment.0.1 += charge.w * tensor.0.1
        octupoleMoment.0.2 += charge.w * tensor.0.2
        // octupoleMoment.1.0 += charge.w * tensor.1.0
        octupoleMoment.1.1 += charge.w * tensor.1.1
        octupoleMoment.1.2 += charge.w * tensor.1.2
        // octupoleMoment.2.0 += charge.w * tensor.2.0
        // octupoleMoment.2.1 += charge.w * tensor.2.1
        octupoleMoment.2.2 += charge.w * tensor.2.2
      }
      
      // Elimination of plane Z, row X.
      octupoleMoment.0.0[2] *= 2
      octupoleMoment.0.1[2] *= 2
      octupoleMoment.0.2 *= 2
      
      // Elimination of plane Z, row Y.
      octupoleMoment.0.1[2] *= 3 / 2
      octupoleMoment.1.1[2] *= 2
      octupoleMoment.1.2 *= 2
      
      // Elimination of plane Y, row X.
      octupoleMoment.0.0[1] *= 2
      octupoleMoment.0.1[1] *= 2
      octupoleMoment.0.1[2] *= 4 / 3
      
      // Partial elimination of row Z in each plane.
      octupoleMoment.0.0[2] *= 3 / 2
      octupoleMoment.0.1[2] *= 5 / 4
      octupoleMoment.0.1[2] *= 6 / 5
      octupoleMoment.1.1[2] *= 3 / 2
      octupoleMoment.0.2 *= 3 / 2
      octupoleMoment.1.2 *= 3 / 2
      
      // Redundantly compute the first element of row Y for simplicity.
      
      let T0_rHat = SIMD3(
        (octupoleMoment.0.0 * rHat).sum(),
        (octupoleMoment.0.1 * rHat).sum(),
        octupoleMoment.0.2 * rHat.z)
      let T1_rHat = SIMD3(
        0,
        (octupoleMoment.1.1 * rHat).sum(),
        octupoleMoment.1.2 * rHat.z)
      let T2_rHat = SIMD3(
        0,
        0,
        octupoleMoment.2.2 * rHat.z)
      let rHat_T_rHat = SIMD3(
        (rHat * T0_rHat).sum(),
        (rHat * T1_rHat).sum(),
        (rHat * T2_rHat).sum())
      let tensorProduct = 1 / 6 * (rHat * rHat_T_rHat).sum()
      potential += tensorProduct / (rNorm * rNorm * rNorm * rNorm)
    }
    estimates.append(potential)
    
    // Evaluate the potential directly.
    potential = .zero
    for charge in pointCharges {
      let position = unsafeBitCast(charge, to: SIMD3<Float>.self)
      let delta = samplePosition - position
      let deltaNorm = (delta * delta).sum().squareRoot()
      potential += charge.w / deltaNorm
    }
    estimates.append(potential)
    
    // Assert that each successive expansion has better quality.
    let actual = estimates.last!
    XCTAssertEqual(estimates[0], actual, accuracy: 0.83)
    XCTAssertEqual(estimates[1], actual, accuracy: 0.13)
    XCTAssertEqual(estimates[2], actual, accuracy: 0.06)
    XCTAssertEqual(estimates[3], actual, accuracy: 0.05)
    XCTAssertEqual(estimates[4], actual, accuracy: 0.01)
    
    XCTAssertNotEqual(estimates[0], actual, accuracy: 0.82)
    XCTAssertNotEqual(estimates[1], actual, accuracy: 0.12)
    XCTAssertNotEqual(estimates[2], actual, accuracy: 0.05)
    XCTAssertNotEqual(estimates[3], actual, accuracy: 0.04)
    XCTAssertNotEqual(estimates[4], actual, accuracy: 0.00)
  }
  
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
  
  // Analyze the case where a cell overlaps itself, and the explicit integral
  // for Hartree potential evaluates to infinity.
  //
  // The 1D case was calculated analytically, and it supposedly diverges.
  // - In the test below, the average value does eventually increase to
  //   infinity. The value scales with the logarithm of the number of grid
  //   points, but log(infinity) is still infinity.
  // The 2D case was calculated numerically: 2.9731896
  // The 3D case was calculated analytically: 2.3800774
  // Source: https://doi.org/10.1103/PhysRevB.50.11355
  //
  // v(r) = Σ ρ(r') g(r, r')
  // ijk ≠ i'j'k' | g(r, r') = h^3 / |r - r'|
  // ijk = i'j'k' | g(r, r') = -h^2 (π / 2 + 3 ln ((√3 - 1) / (√3 + 1)))
  func testSelfRepulsionIntegral1D() throws {
    // Try integrating the 1D and 2D integrals numerically.
    // 1) Exclude the singularities.
    // 2) Find the limit as the singularity becomes infinitesimally small.
    // 3) Take the average over the entire volume.
    let h: Float = 0.004
    let cellCount: Int = 250
    
    // Create an array that represents the charge density.
    //
    // Change the default value to check that the potential solver is correctly
    // responding to the input.
    let chargeDistribution = [Float](repeating: 1, count: cellCount)
    
    // Solve for the potential distribution with an integral.
    var potentialDistribution: [Float] = []
    for cellID in 0..<cellCount {
      let x = (Float(cellID) + 0.5) * h
      
      var accumulator: Double = .zero
      for otherCellID in 0..<cellCount where cellID != otherCellID {
        let otherX = (Float(otherCellID) + 0.5) * h
        
        // Integrate over 1D space.
        let ρ = chargeDistribution[otherCellID]
        let g = h / (x - otherX).magnitude
        accumulator += Double(ρ * g)
      }
      
      let v = Float(accumulator)
      potentialDistribution.append(v)
    }
    
    // Visualize the potential.
    do {
      let potentialLeft = potentialDistribution[0]
      let potentialLeft2 = potentialDistribution[1]
      let potentialCenter = potentialDistribution[cellCount / 2]
      let potentialRight2 = potentialDistribution[cellCount - 2]
      let potentialRight = potentialDistribution[cellCount - 1]
      XCTAssertEqual(potentialLeft, 6.0966754, accuracy: 1e-3)
      XCTAssertEqual(potentialLeft2, 7.092659, accuracy: 1e-3)
      XCTAssertEqual(potentialCenter, 10.811043, accuracy: 1e-3)
      XCTAssertEqual(potentialRight2, 7.0926757, accuracy: 1e-3)
      XCTAssertEqual(potentialRight, 6.0966783, accuracy: 1e-3)
    }
    
    // Report the average of the potential over the entire domain.
    do {
      var accumulator: Double = .zero
      for cellID in 0..<cellCount {
        let v = potentialDistribution[cellID]
        let drTerm = h
        accumulator += Double(v * drTerm)
      }
      
      let average = Float(accumulator)
      XCTAssertEqual(average, 10.201351, accuracy: 1e-3)
    }
  }
  
  // TODO: Clean up this test case, turn into a proper unit test.
  func testSelfRepulsionIntegral2D() throws {
    let h: Float = 1.0 / 10
    let gridSize: Int = 10
    
    // Create an array that represents the charge density.
    let chargeGrid = [Float](repeating: 1, count: gridSize * gridSize)
    
    // Solve for the potential distribution with an integral.
    var potentialGrid: [Float] = []
    for indexY in 0..<gridSize {
      for indexX in 0..<gridSize {
        let x = (Float(indexX) + 0.5) * h
        let y = (Float(indexY) + 0.5) * h
        let cellID = indexY * gridSize + indexX
        
        var accumulator: Double = .zero
        for otherIndexY in 0..<gridSize {
          for otherIndexX in 0..<gridSize {
            let otherX = (Float(otherIndexX) + 0.5) * h
            let otherY = (Float(otherIndexY) + 0.5) * h
            let otherCellID = otherIndexY * gridSize + otherIndexX
            
            // Skip the singularity.
            if cellID == otherCellID {
              continue
            }
            
            // Find the distance between the two points.
            let r = SIMD2(x, y)
            let rPrime = SIMD2(otherX, otherY)
            let rDelta = r - rPrime
            let distance = (rDelta * rDelta).sum().squareRoot()
            
            // Integrate over 2D space.
            let ρ = chargeGrid[otherCellID]
            let g = (h * h) / distance
            accumulator += Double(ρ * g)
          }
        }
        
        let v = Float(accumulator)
        potentialGrid.append(v)
      }
    }
    
    // Visualize the potential.
    do {
      let size = gridSize
      
      let potential00 = potentialGrid[0 * size + 0]
      let potential01 = potentialGrid[0 * size + (size / 2)]
      let potential02 = potentialGrid[0 * size + (size - 1)]
      print(potential00, potential01, potential02)
      
      let potential10 = potentialGrid[(size / 2) * size + 0]
      let potential11 = potentialGrid[(size / 2) * size + (size / 2)]
      let potential12 = potentialGrid[(size / 2) * size + (size - 1)]
      print(potential10, potential11, potential12)
      
      let potential20 = potentialGrid[(size - 1) * size + 0]
      let potential21 = potentialGrid[(size - 1) * size + (size / 2)]
      let potential22 = potentialGrid[(size - 1) * size + (size - 1)]
      print(potential20, potential21, potential22)
    }
    
    // Report the average of the potential over the entire domain.
    //
    // This sequence converges:
    // size =   5 | 2.2688484 |
    // size =  10 | 2.6065936 | +0.3377452
    // size =  20 | 2.785188 | +0.1785944
    // size =  40 | 2.8777378 | +0.0925498
    // size =  80 | 2.9250371 | +0.0472993
    // size = 160 | 2.9489965 | +0.0239594
    // size = 320 | 2.961067 | +0.0120705
    // size = 640 | 2.9671283 | +0.0060613
    // size = inf | 2.9731896 |
    do {
      var accumulator: Double = .zero
      for cellID in 0..<gridSize * gridSize {
        let v = potentialGrid[cellID]
        let drTerm = h * h
        accumulator += Double(v * drTerm)
      }
      
      let average = Float(accumulator)
      print()
      print(average)
    }
  }
}
