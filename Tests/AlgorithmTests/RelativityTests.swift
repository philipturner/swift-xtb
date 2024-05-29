import XCTest

final class RelativityTests: XCTestCase {
  // Test the recursive definition of kinetic energy.
  // Source: https://doi.org/10.1088/1361-6404/ac0ecc
  //
  // Key formulas:
  // γ = √(1 + v^2 / c^2)
  // K_Sch = p^2 / (2m)
  // K_GP = p^2 / ((γ + 1)m)
  // K_P = p^2 / (2m + K_P / c^2))
  func testKineticEnergyFormula() throws {
    func testQuasiRelativisticEnergy(v: Double) {
      // Momentum and velocity are the same in atomic units, as mass is 1.
      let p: Double = v
      let c: Double = 137.035999084
      
      // The direct formulas for kinetic energy.
      let gamma = Double(1 + p * p / (c * c)).squareRoot()
      let kineticEnergySch: Double = p * p / 2
      let kineticEnergyGP = p * p / (gamma + 1)
      
      // The recursive formula for kinetic energy.
      var kineticEnergyP = kineticEnergySch
      for _ in 0..<15 {
        kineticEnergyP = p * p / (2 + kineticEnergyP / (c * c))
      }
      XCTAssertEqual(kineticEnergyGP, kineticEnergyP, accuracy: 0.001)
    }
    
    testQuasiRelativisticEnergy(v: 1)
    testQuasiRelativisticEnergy(v: 14)
    testQuasiRelativisticEnergy(v: 32)
    testQuasiRelativisticEnergy(v: 79)
    testQuasiRelativisticEnergy(v: 137)
    testQuasiRelativisticEnergy(v: 172)
  }
  
  // Solve the radial wavefunction for hydrogen-like atoms. Compare the results
  // to the paper where BigDFT did quasi-relativistic quantum mechanics.
  // - Does shifting the kinetic energy change the amplitude of the
  //   wavefunction anywhere? The wavefunction must change so that other
  //   electrons can see the effects.
  // - Log the amplitude of the wavefunction at certain points, with and
  //   without the relativistic correction. Also, include the energies with the
  //   wavefunctions.
  //
  // Data for energies of different atoms, with different approximations. These
  // pertain to the 1s orbital of the one-electron system.
  // Source: https://doi.org/10.1016/j.comptc.2020.112711
  // Source: https://doi.org/10.1063/1.1818681
  //
  // |  Z  |  E (Sch)  |   E (RK)  |  E (DKH1) |  E (DEQ)  |  E (PGPP) |
  // | --- | --------- | --------- | --------- | --------- | --------- |
  // |   1 |    -0.500 |    -0.500 |    -0.500 |       n/a |    -0.500 |
  // |   8 |   -32.000 |   -32.128 |   -32.031 |       n/a |   -32.031 |
  // |  20 |  -200.000 |  -204.807 |  -201.341 |  -201.077 |  -201.091 |
  // |  40 |  -800.000 |  -878.142 |  -823.894 |  -817.807 |  -817.502 |
  // |  60 | -1800.000 | -2252.407 | -1934.203 | -1895.682 | -1890.845 |
  // |  80 | -3200.000 | -5317.788 | -3686.447 | -3532.192 | -3498.570 |
  //
  // Insights: the results of DEQ and PGPP agree remarkably well. DEQ is the
  // "exact" version of DKH, and its results were often very close to DKH14 in
  // the paper. However, this test only considers the one-electron case, and
  // neglects spin-orbit coupling. PGPP could perform significantly worse for
  // many-electron systems.
  //
  // Side-by-side comparison of wavefunction values and energies for Au(79):
  //
  //                     | Nonrelativistic | Relativistic  |
  // ------------------- | --------------- | ------------- |
  // Energy              | -3120.803       | -3403.7566    |
  // <r>                 | 0.012657344     | 0.0116049675  |
  // ------------------- | --------------- | ------------- |
  // Ψ(0.000025 * 0.5)   | 392.45526       | 446.65286     |
  // Ψ(0.000025 * 5.5)   | 355.47208       | 400.93982     |
  // Ψ(0.000025 * 10.5)  | 322.02545       | 359.9737      |
  // Ψ(0.000025 * 50.5)  | 146.12407       | 152.04948     |
  // Ψ(0.000025 * 100.5) | 54.42798        | 51.78401      |
  // Ψ(0.000025 * 200.5) | 7.5521126       | 6.0070696     |
  // Ψ(0.000025 * 300.5) | 1.0479369       | 0.6968856     |
  // Ψ(0.000025 * 400.5) | 0.14541396      | 0.08085826    |
  // Ψ(0.000025 * 500.5) | 0.020177867     | 0.009385774   |
  // Ψ(0.000025 * 600.5) | 0.002799893     | 0.0010906607  |
  // Ψ(0.000025 * 700.5) | 0.000388516     | 0.00012704461 |
  // Ψ(0.000025 * 800.5) | 5.3894604e-05   | 1.4866069e-05 |
  // Ψ(0.000025 * 900.5) | 7.3126744e-06   | 1.7254465e-06 |
  // Ψ(0.000025 * 999.5) | 3.9033416e-08   | 8.479288e-09  |
  //
  // The radius only contracts by 1.09x. This is less than the Wikipedia value
  // of 1.22x. However, the number on Wikipedia (and found on another lecture
  // PDF), is a gross approximation based on the Lorentz factor. It does not
  // enter the formula into an SCF calculation and report the result.
  func testHydrogenicAtom() throws {
    // Settings: cellCount = 1000, SCF iterations = 1000, nonrelativistic
    // h = 0.02, Z = 1 -> -0.50005084
    // h = 0.0025, Z = 8 -> -32.00325
    // h = 0.001, Z = 20 -> -200.01988
    // h = 0.0005, Z = 40 -> -800.0795
    // h = 0.000333, Z = 60 -> -1800.1755
    // h = 0.00025, Z = 80 -> -3200.318
    let h: Float = 0.001
    let cellCount: Int = 250
    let scfIterationCount: Int = 200
    let Z: Int = 79
    
    // Applies the Laplacian operator to the wavefunction.
    func laplacian(Ψ: [Float]) -> [Float] {
      var output: [Float] = []
      for cellID in 0..<cellCount {
        let r = (Float(cellID) + 0.5) * h
        
        // Gather data points for the finite difference.
        var left: Float = .zero
        var center: Float = .zero
        var right: Float = .zero
        if Ψ.indices.contains(cellID - 1) {
          left = Ψ[cellID - 1]
        } else {
          // The wavefunction is spherical about the origin.
          left = Ψ[cellID]
        }
        if Ψ.indices.contains(cellID) {
          center = Ψ[cellID]
        }
        if Ψ.indices.contains(cellID + 1) {
          right = Ψ[cellID + 1]
        }
        
        // Evaluate the Laplacian in spherical coordinates.
        // Source: https://en.wikipedia.org/wiki/Laplace_operator#Three_dimensions
        
        // Apply the derivative the first time.
        var leftDerivative = (center - left) / h
        var rightDerivative = (right - center) / h
        
        // Scale by r^2.
        leftDerivative *= (r - 0.5 * h) * (r - 0.5 * h)
        rightDerivative *= (r + 0.5 * h) * (r + 0.5 * h)
        
        // Apply the derivative the second time.
        var secondDerivative = (rightDerivative - leftDerivative) / h
        
        // Undo the scaling by r^2.
        secondDerivative /= (r * r)
        output.append(secondDerivative)
      }
      return output
    }
    
    // Applies the Hamiltonian operator to the wavefunction.
    func hamiltonian(Ψ: [Float], γ: Float) -> [Float] {
      let laplacianΨ = laplacian(Ψ: Ψ)
      
      var output: [Float] = []
      for cellID in 0..<cellCount {
        let r = (Float(cellID) + 0.5) * h
        let amplitude = Ψ[cellID]
        
        // Evaluate the energy terms.
        let LΨ = laplacianΨ[cellID]
        let TΨ = -1 / (1 + γ) * LΨ
        let VΨ = -Float(Z) / r * amplitude
        
        // Return the energy times the wavefunction.
        let HΨ = TΨ + VΨ
        output.append(HΨ)
      }
      return output
    }
    
    // Finds the relativistic gamma factor.
    func createGamma(Ψ: [Float]) -> Float {
      let laplacianΨ = laplacian(Ψ: wavefunction)
      let ΨLΨ = integral(wavefunction, laplacianΨ)
      let p2 = -ΨLΨ
      
      let c: Double = 137.035999084
      var γ = 1 + Double(p2) / (c * c)
      γ.formSquareRoot()
      return Float(γ)
    }
    
    // Evaluates an integral in spherical coordinates.
    func integral(_ lhs: [Float], _ rhs: [Float]) -> Float {
      guard lhs.count == cellCount,
            rhs.count == cellCount else {
        fatalError("Operands must have same size as domain.")
      }
      
      // Accumulate in FP64, then convert the sum to FP32.
      var accumulator: Double = .zero
      for cellID in 0..<cellCount {
        let r = (Float(cellID) + 0.5) * h
        let product = lhs[cellID] * rhs[cellID]
        let drTerm = 4 * Float.pi * (r * r) * h
        accumulator += Double(product * drTerm)
      }
      return Float(accumulator)
    }
    
    // Shifts the LHS by a factor times the RHS.
    func shift(
      _ lhs: [Float], _ factor: Float, _ rhs: [Float]
    ) -> [Float] {
      guard lhs.count == cellCount,
            rhs.count == cellCount else {
        fatalError("Operands must have same size as domain.")
      }
      
      var output: [Float] = []
      for cellID in 0..<cellCount {
        let lhsValue = lhs[cellID]
        let rhsValue = rhs[cellID]
        let outValue = lhsValue + factor * rhsValue
        output.append(outValue)
      }
      return output
    }
    
    // Hydrogen wave function: R(r) = 2 e^{-r}
    //
    // Solve the wavefunction on a radial grid for simplicity.
    var wavefunction: [Float] = []
    for cellID in 0..<cellCount {
      let r = (Float(cellID) + 0.5) * h
      let R = 2 * Float.exp(-Float(Z) * r)
      
      var normalizationFactor = 1 / (4 * Float.pi)
      normalizationFactor.formSquareRoot()
      normalizationFactor *= Float.pow(Float(Z), 1.5)
      wavefunction.append(normalizationFactor * R)
    }
    
    // Check the normalization factor.
    do {
      let normalizationFactor = integral(wavefunction, wavefunction)
      XCTAssertEqual(normalizationFactor, 1, accuracy: 1e-3)
    }
    
    // Report the energy before doing the SCF iterations.
    do {
      let γ = createGamma(Ψ: wavefunction)
      let hamiltonianΨ = hamiltonian(Ψ: wavefunction, γ: γ)
      let ΨHΨ = integral(wavefunction, hamiltonianΨ)
      let ΨΨ  = integral(wavefunction, wavefunction)
      let rayleighQuotient = ΨHΨ / ΨΨ
      XCTAssertEqual(rayleighQuotient, -3349.0366, accuracy: 1e-1)
    }
    
    // Search for the lowest eigenpair.
    for _ in 0..<scfIterationCount {
      // Cache the hamiltonian times the wavefunction.
      let γ = createGamma(Ψ: wavefunction)
      var hamiltonianΨ = hamiltonian(Ψ: wavefunction, γ: γ)
      
      for _ in 0..<5 {
        // Find the Rayleigh quotient.
        let ΨHΨ = integral(wavefunction, hamiltonianΨ)
        let ΨΨ  = integral(wavefunction, wavefunction)
        let rayleighQuotient = ΨHΨ / ΨΨ
        let residual = shift(hamiltonianΨ, -rayleighQuotient, wavefunction)
        
        // Evaluate the integrals necessary to construct the timestep.
        let hamiltonianR = hamiltonian(Ψ: residual, γ: γ)
        let rr = integral(residual, residual)
        let Ψr = integral(wavefunction, residual)
        let rHr = integral(residual, hamiltonianR)
        let ΨHr = integral(wavefunction, hamiltonianR)
        
        // Solve a quadratic equation.
        let m = [rr, Ψr, rHr, ΨHr, ΨHΨ, ΨΨ]
        let a = m[0] * m[3] - m[2] * m[1]
        let b = m[5] * m[2] - m[4] * m[0]
        let c = m[4] * m[1] - m[3] * m[5]
        let den = b + (b * b - 4 * a * c).squareRoot()
        
        // Choose the timestep.
        var λ: Float = .zero
        if den.magnitude > 1e-15 {
          λ = 2 * c / den
        }
        
        // Update the wavefunction.
        wavefunction = shift(wavefunction, λ, residual)
        hamiltonianΨ = shift(hamiltonianΨ, λ, hamiltonianR)
      }
      
      // Normalize the wavefunction.
      let norm = integral(wavefunction, wavefunction)
      let scaleFactor = 1 / norm.squareRoot()
      for cellID in 0..<cellCount {
        wavefunction[cellID] *= scaleFactor
      }
    }
    
    // Utility function to check that the calculation is converging.
    func assertProgress(start: Float, current: Float, end: Float) {
      let progress = (current - start) / (end - start)
      XCTAssertGreaterThan(
        progress, 0.50, "Less than 50% way through convergence.")
      XCTAssertLessThan(
        progress, 1.20, "More than 120% way through convergence.")
    }
    
    // Report the energy after doing the SCF iterations.
    do {
      let γ = createGamma(Ψ: wavefunction)
      let hamiltonianΨ = hamiltonian(Ψ: wavefunction, γ: γ)
      let ΨHΨ = integral(wavefunction, hamiltonianΨ)
      let ΨΨ  = integral(wavefunction, wavefunction)
      let rayleighQuotient = ΨHΨ / ΨΨ
      assertProgress(
        start: -3349.0366, current: rayleighQuotient, end: -3411.3223)
    }
    
    // Report the expectation value for radius.
    do {
      // Multiply the radius with the wavefunction.
      var rΨ: [Float] = []
      for cellID in 0..<cellCount {
        let r = (Float(cellID) + 0.5) * h
        let amplitude = wavefunction[cellID]
        rΨ.append(r * amplitude)
      }
      
      // Find the expectation value by taking an integral. As shown in previous
      // experiments, the radius needs to be multiplied by 2/3 to get the
      // "atomic radius".
      let ΨrΨ = integral(wavefunction, rΨ)
      let atomicRadius = ΨrΨ * 2.0 / 3
      assertProgress(
        start: 0.012658219, current: atomicRadius, end: 0.011579066)
    }
    
    // Report the wavefunction values.
    do {
      var indices: [Int] = []
      indices.append(0)
      indices.append(1)
      indices.append(2)
      indices.append(10)
      indices.append(20)
      indices.append(50)
      indices.append(100)
      indices.append(150)
      indices.append(200)
      indices.append(249)
      
      var startAmplitudes: [Float] = []
      startAmplitudes.append(380.8124)
      startAmplitudes.append(351.88586)
      startAmplitudes.append(325.15656)
      startAmplitudes.append(172.82971)
      startAmplitudes.append(78.43787)
      startAmplitudes.append(7.3324285)
      startAmplitudes.append(0.14118369)
      startAmplitudes.append(0.0027184465)
      startAmplitudes.append(5.2342937e-05)
      startAmplitudes.append(1.0906967e-06)
      XCTAssertEqual(indices.count, startAmplitudes.count)
      
      var endAmplitudes: [Float] = []
      endAmplitudes.append(435.51288)
      endAmplitudes.append(398.7813)
      endAmplitudes.append(365.48648)
      endAmplitudes.append(182.82532)
      endAmplitudes.append(77.071526)
      endAmplitudes.append(5.78571)
      endAmplitudes.append(0.07740794)
      endAmplitudes.append(0.0010362118)
      endAmplitudes.append(1.3870742e-05)
      endAmplitudes.append(3.0690966e-08)
      XCTAssertEqual(indices.count, endAmplitudes.count)
      
      for indexID in indices.indices {
        let startAmplitude = startAmplitudes[indexID]
        let index = indices[indexID]
        let endAmplitude = endAmplitudes[indexID]
        
        let amplitude = wavefunction[index]
        assertProgress(
          start: startAmplitude, current: amplitude, end: endAmplitude)
      }
    }
  }
}
