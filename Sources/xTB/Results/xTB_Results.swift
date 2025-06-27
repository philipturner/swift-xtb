//
//  xTB_Results.swift
//  swift-xtb
//
//  Created by Philip Turner on 5/30/24.
//

import C_xTB

class xTB_Results {
  unowned var calculator: xTB_Calculator!
  
  var tResults: xtb_TResults!
  
  var energy: Double?
  
  var forces: [SIMD3<Float>]?
  
  var charges: [Float]?
  
  var bondOrders: [Float]?
  
  var orbitalEigenvalues: [Float]?
  
  var orbitalOccupations: [Float]?
  
  var orbitalCoefficients: [Float]?
  
  init() {
    self.tResults = xTB_Results.createObject()
  }
  
  deinit {
    xtb_delResults(&tResults)
  }
  
  /// Create the reference-counted object from the C API.
  static func createObject() -> xtb_TResults {
    let res = xtb_newResults()
    guard let res else {
      fatalError("Could not create new xTB_Results.")
    }
    return res
  }
  
  private typealias DoubleArrayFunction = @convention(c) (
    xtb_TEnvironment?,
    xtb_TResults?,
    UnsafeMutablePointer<Double>?
  ) -> Void
  
  private func getDoubleArray(
    symbol: DoubleArrayFunction,
    size: Int
  ) -> [Double] {
    var output = [Double](repeating: .zero, count: size)
    symbol(xTB_Environment.tEnvironment, tResults, &output)
    return output
  }
}

// MARK: - Energy

extension xTB_Results {
  func getEnergy() -> Double {
    var energy: Double = .zero
    xtb_getEnergy(
      xTB_Environment.tEnvironment, tResults, &energy)
    
    // Convert energy into nanomechanical units.
    return energy * xTB_ZJPerHartree
  }
}

// MARK: - Molecule

extension xTB_Results {
  func getForces() -> [SIMD3<Float>] {
    let atomCount = calculator.molecule.atomicNumbers.count
    let gradient64 = getDoubleArray(
      symbol: xtb_getGradient,
      size: atomCount * 3)
    
    // Convert forces into nanomechanical units and flip their sign.
    return convertGradientToForces(gradient64)
  }
  
  func getCharges() -> [Float] {
    let atomCount = calculator.molecule.atomicNumbers.count
    let charges64 = getDoubleArray(
      symbol: xtb_getCharges,
      size: atomCount)
    return charges64.map(Float.init)
  }
  
  func getBondOrders() -> [Float] {
    let atomCount = calculator.molecule.atomicNumbers.count
    let bondOrders64 = getDoubleArray(
      symbol: xtb_getBondOrders,
      size: atomCount * atomCount)
    return bondOrders64.map(Float.init)
  }
}

// MARK: - Orbitals

extension xTB_Results {
  func checkOrbitalCount() {
    var orbitalCount: Int32 = .max
    xtb_getNao(
      xTB_Environment.tEnvironment, tResults, &orbitalCount)
    guard calculator.orbitals.count == Int(orbitalCount) else {
      fatalError("Orbital count did not match expectations.")
    }
  }
  
  func getOrbitalEigenvalues() -> [Float] {
    switch calculator.hamiltonian {
    case .forceField:
      return []
    case .tightBinding:
      let orbitalCount = calculator.orbitals.count
      let orbitalEigenvalues64 = getDoubleArray(
        symbol: xtb_getOrbitalEigenvalues,
        size: orbitalCount)
      
      // Convert energy into nanomechanical units.
      return orbitalEigenvalues64.map {
        Float($0) * Float(xTB_ZJPerHartree)
      }
    }
  }
  
  func getOrbitalOccupations() -> [Float] {
    switch calculator.hamiltonian {
    case .forceField:
      return []
    case .tightBinding:
      let orbitalCount = calculator.orbitals.count
      let orbitalOccupations64 = getDoubleArray(
        symbol: xtb_getOrbitalOccupations,
        size: orbitalCount)
      return orbitalOccupations64.map(Float.init)
    }
  }
  
  func getOrbitalCoefficients() -> [Float] {
    switch calculator.hamiltonian {
    case .forceField:
      return []
    case .tightBinding:
      let orbitalCount = calculator.orbitals.count
      let orbitalCoefficients64 = getDoubleArray(
        symbol: xtb_getOrbitalCoefficients,
        size: orbitalCount * orbitalCount)
      return orbitalCoefficients64.map(Float.init)
    }
  }
}
