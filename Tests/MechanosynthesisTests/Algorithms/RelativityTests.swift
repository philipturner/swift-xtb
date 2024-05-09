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
  func testHydrogenicAtom() throws {
    let h: Float = 0.1
    let cellCount: Int = 100
    let Z: Int = 1
    
    // Hydrogen wave function: R(r) = 2 e^{-r}
    //
    // Solve the wavefunction on a radial grid for simplicity. Test the
    // expectation value for normalization factor and energy.
    // - Ensure the naive version integrates to 4 pi.
    var wavefunction: [Float] = []
    for cellID in 0..<cellCount {
      let r = (Float(cellID) + 0.5) * h
      let radialTerm = 2 * Float.exp(-r)
      
      var normalizationFactor = 1 / (4 * Float.pi)
      normalizationFactor.formSquareRoot()
      wavefunction.append(normalizationFactor * radialTerm)
    }
    
    // Check the normalization factor.
    do {
      var normalizationFactorSum: Double = .zero
      for cellID in 0..<cellCount {
        let amplitude = wavefunction[cellID]
        let valueTerm = amplitude * amplitude
        
        let r = (Float(cellID) + 0.5) * h
        let drTerm = 4 * Float.pi * (r * r) * h
        normalizationFactorSum += Double(valueTerm * drTerm)
      }
      let normalizationFactor = Float(normalizationFactorSum)
      XCTAssertEqual(normalizationFactor, 1, accuracy: 1e-5)
    }
    
    // Query the energy of the ansatz.
    do {
      var kineticEnergySum: Double = .zero
      var potentialEnergySum: Double = .zero
      for cellID in 0..<cellCount {
        let r = (Float(cellID) + 0.5) * h
        
        func createLaplacian() -> Float {
          // Gather data points for the finite difference.
          var left: Float = .zero
          var center: Float = .zero
          var right: Float = .zero
          if wavefunction.indices.contains(cellID - 1) {
            left = wavefunction[cellID - 1]
          } else {
            // The wavefunction is spherical about the origin.
            left = wavefunction[cellID]
          }
          if wavefunction.indices.contains(cellID) {
            center = wavefunction[cellID]
          }
          if wavefunction.indices.contains(cellID + 1) {
            right = wavefunction[cellID + 1]
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
          return secondDerivative
        }
        
        // Fetch the value of the wavefunction.
        let amplitude = wavefunction[cellID]
        
        // Evaluate the energy terms.
        let laplacian = createLaplacian()
        var kineticEnergyTerm = -0.5 * createLaplacian()
        var potentialEnergyTerm = -Float(Z) / r * amplitude
        kineticEnergyTerm *= amplitude
        potentialEnergyTerm *= amplitude
        
        // Sum the terms into the integral.
        let drTerm = 4 * Float.pi * (r * r) * h
        kineticEnergySum += Double(kineticEnergyTerm * drTerm)
        potentialEnergySum += Double(potentialEnergyTerm * drTerm)
        print(String(format: "%.3f", r), laplacian, kineticEnergyTerm, potentialEnergyTerm)
      }
      
      // Log the values returned by the integral.
      let kineticEnergy = Float(kineticEnergySum)
      let potentialEnergy = Float(potentialEnergySum)
      print()
      print(kineticEnergy)
      print(potentialEnergy)
    }
  }
}
