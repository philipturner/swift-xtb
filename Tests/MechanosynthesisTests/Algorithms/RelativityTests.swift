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
}
