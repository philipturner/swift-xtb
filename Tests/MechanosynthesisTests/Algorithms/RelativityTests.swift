import XCTest
import Mechanosynthesis
import Numerics

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
  // |  Z  |  E (Sch)  |   E (RK)  |  E (DKH1) |  E (DEQ)  | This Test |
  // | --- | --------- | --------- | --------- | --------- | --------- |
  // |   1 |    -0.500 |    -0.500 |    -0.500 |       n/a |
  // |   8 |   -32.000 |   -32.128 |   -32.031 |       n/a |
  // |  20 |  -200.000 |  -204.807 |  -201.341 |  -201.077 |
  // |  40 |  -800.000 |  -878.142 |  -823.894 |  -817.807 |
  // |  60 | -1800.000 | -2252.407 | -1934.203 | -1895.682 |
  // |  80 | -3200.000 | -5317.788 | -3686.447 | -3532.192 |
  //
  // Side-by-side comparison of wavefunction values and energies for Au(79):
  //
  //                    | Nonrelativistic | Relativistic |
  // ------------------ | --------------- | ------------ |
  // Energy             |
  // Ψ(0.000025 * 0.5)  |
  // Ψ(0.000025 * 5.5)  |
  // Ψ(0.000025 * 10.5) |
  func testHydrogenicAtom() throws {
    // Nonrelativistic results (for debugging):
    //
    // Settings: cellCount = 1000, SCF iterations = 1000
    // h = 0.02, Z = 1 -> -0.50005084
    // h = 0.0025, Z = 8 -> -32.00325
    // h = 0.001, Z = 20 -> -200.01988
    // h = 0.0005, Z = 40 -> -800.0795
    // h = 0.000333, Z = 60 -> -1800.1755
    // h = 0.00025, Z = 80 -> -3200.318
    let h: Float = 0.02
    let cellCount: Int = 1000
    let Z: Int = 1
    
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
      var hamiltonianΨ = hamiltonian(Ψ: wavefunction, γ: 1)
      let ΨHΨ = integral(wavefunction, hamiltonianΨ)
      let ΨΨ  = integral(wavefunction, wavefunction)
      let rayleighQuotient = ΨHΨ / ΨΨ
      print()
      print("initial energy:", rayleighQuotient)
    }
    
    // Search for the lowest eigenpair.
    print()
    for _ in 0..<2 {
      // Cache the hamiltonian times the wavefunction.
      let γ = createGamma(Ψ: wavefunction)
      var hamiltonianΨ = hamiltonian(Ψ: wavefunction, γ: 1)
      
      for iterationID in 0..<5 {
        // Find the Rayleigh quotient.
        let ΨHΨ = integral(wavefunction, hamiltonianΨ)
        let ΨΨ  = integral(wavefunction, wavefunction)
        let rayleighQuotient = ΨHΨ / ΨΨ
        let residual = shift(hamiltonianΨ, -rayleighQuotient, wavefunction)
        
        // Evaluate the integrals necessary to construct the timestep.
        let hamiltonianR = hamiltonian(Ψ: residual, γ: 1)
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
        if iterationID == 4 {
          print("normres:", rr.squareRoot(), "| λ:", λ, "| γ:", γ)
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
    
    // Report the energy after doing the SCF iterations.
    do {
      var hamiltonianΨ = hamiltonian(Ψ: wavefunction, γ: 1)
      let ΨHΨ = integral(wavefunction, hamiltonianΨ)
      let ΨΨ  = integral(wavefunction, wavefunction)
      let rayleighQuotient = ΨHΨ / ΨΨ
      print()
      print("final energy:", rayleighQuotient)
    }
  }
}
