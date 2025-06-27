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
    symbol(
      xTB_Environment.tEnvironment,
      calculator.results.tResults, // guarantee source of truth
      &output)
    return output
  }
}

// MARK: - Energy

extension xTB_Results {
  func getEnergy() {
    print("xtb_getEnergy")
    var energy: Double = .zero
    xtb_getEnergy(
      xTB_Environment.tEnvironment,
      calculator.results.tResults, // guarantee source of truth
      &energy)
    
    // Convert energy into nanomechanical units.
    self.energy = energy * xTB_ZJPerHartree
  }
}

// MARK: - Molecule

extension xTB_Results {
  func getForces() {
    print("xtb_getForces")
    let atomCount = calculator.molecule.atomicNumbers.count
    let gradient64 = getDoubleArray(
      symbol: xtb_getGradient,
      size: atomCount * 3)
    forces = convertGradientToForces(gradient64)
  }
  
  func getCharges() {
    print("xtb_getCharges")
    let atomCount = calculator.molecule.atomicNumbers.count
    let charges64 = getDoubleArray(
      symbol: xtb_getCharges,
      size: atomCount)
    charges = charges64.map(Float.init)
  }
  
  func getBondOrders() {
    print("xtb_getBondOrders")
    let atomCount = calculator.molecule.atomicNumbers.count
    let bondOrders64 = getDoubleArray(
      symbol: xtb_getBondOrders,
      size: atomCount * atomCount)
    bondOrders = bondOrders64.map(Float.init)
  }
}

// MARK: - Orbitals

extension xTB_Results {
  func checkOrbitalCount() {
    print("xtb_getNao")
    var orbitalCount: Int32 = .max
    xtb_getNao(
      xTB_Environment.tEnvironment,
      calculator.results.tResults, // guarantee source of truth
      &orbitalCount)
    guard calculator.orbitals.count == Int(orbitalCount) else {
      fatalError("Orbital count did not match expectations.")
    }
  }
  
  func getOrbitalEigenvalues() {
    print("xtb_getOrbitalEigenvalues")
    let orbitalCount = calculator.orbitals.count
    let orbitalEigenvalues64 = getDoubleArray(
      symbol: xtb_getOrbitalEigenvalues,
      size: orbitalCount)
    
    // Convert energy into nanomechanical units.
    orbitalEigenvalues = orbitalEigenvalues64.map {
      Float($0) * Float(xTB_ZJPerHartree)
    }
  }
  
  func getOrbitalOccupations() {
    print("xtb_getOrbitalOccupations")
    let orbitalCount = calculator.orbitals.count
    let orbitalOccupations64 = getDoubleArray(
      symbol: xtb_getOrbitalOccupations,
      size: orbitalCount)
    orbitalOccupations = orbitalOccupations64.map(Float.init)
  }
  
  func getOrbitalCoefficients() {
    print("xtb_getOrbitalCoefficients")
    let orbitalCount = calculator.orbitals.count
    let orbitalCoefficients64 = getDoubleArray(
      symbol: xtb_getOrbitalCoefficients,
      size: orbitalCount * orbitalCount)
    orbitalCoefficients = orbitalCoefficients64.map(Float.init)
  }
}
