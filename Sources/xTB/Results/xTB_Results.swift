//
//  xTB_Results.swift
//  
//
//  Created by Philip Turner on 5/30/24.
//

class xTB_Results {
  unowned var calculator: xTB_Calculator!
  
  var _results: xtb_TResults
  
  var energy: Double?
  
  var forces: [SIMD3<Float>]?
  
  var externalChargeForces: [SIMD3<Float>]?
  
  var charges: [Float]?
  
  var bondOrders: [Float]?
  
  var orbitalEigenvalues: [Float]?
  
  var orbitalOccupations: [Float]?
  
  var orbitalCoefficients: [Float]?
  
  init() {
    guard let res = xtb_newResults() else {
      fatalError("Could not create new xTB_Results.")
    }
    _results = res
  }
  
  deinit {
    xtb_delResults(&_results)
  }
}

extension xTB_Results {
  private typealias DoubleArrayFunction = @convention(c) (
    xtb_TEnvironment,
    xtb_TResults,
    UnsafeMutablePointer<Double>?
  ) -> Void
  
  private func getDoubleArray(
    _ closure: DoubleArrayFunction,
    size: Int
  ) -> [Double] {
    var output = [Double](repeating: .zero, count: size)
    closure(
      xTB_Environment._environment,
      calculator.results._results,
      &output)
    return output
  }
}

extension xTB_Results {
  func getEnergy() {
    var energy: Double = .zero
    xtb_getEnergy(
      xTB_Environment._environment,
      calculator.results._results,
      &energy)
    
    // Convert energy into nanomechanical units.
    self.energy = energy * xTB_ZJPerHartree
  }
  
  func getExternalChargeForces() {
    let externalChargeCount = calculator.externalCharges.atomicNumbers.count
    let pointChargeGradient64 = getDoubleArray(
      xtb_getPCGradient, size: externalChargeCount * 3)
    externalChargeForces = convertGradientToForces(pointChargeGradient64)
  }
}

extension xTB_Results {
  func getForces() {
    let atomCount = calculator.molecule.atomicNumbers.count
    let gradient64 = getDoubleArray(
      xtb_getGradient, size: atomCount * 3)
    forces = convertGradientToForces(gradient64)
  }
  
  func getCharges() {
    let atomCount = calculator.molecule.atomicNumbers.count
    let charges64 = getDoubleArray(
      xtb_getCharges, size: atomCount * 3)
    charges = charges64.map(Float.init)
  }
  
  func getBondOrders() {
    let atomCount = calculator.molecule.atomicNumbers.count
    let bondOrders64 = getDoubleArray(
      xtb_getBondOrders, size: atomCount * atomCount)
    bondOrders = bondOrders64.map(Float.init)
  }
}

extension xTB_Results {
  func checkOrbitalCount() {
    var orbitalCount: Int32 = .max
    xtb_getNao(
      xTB_Environment._environment,
      calculator.results._results,
      &orbitalCount)
    guard calculator.orbitals.count == Int(orbitalCount) else {
      fatalError("Orbital count did not match expectations.")
    }
  }
  
  func getOrbitalEigenvalues() {
    let orbitalCount = calculator.orbitals.count
    let orbitalEigenvalues64 = getDoubleArray(
      xtb_getOrbitalEigenvalues, size: orbitalCount)
    
    // Convert energy into nanomechanical units.
    orbitalEigenvalues = orbitalEigenvalues64.map {
      Float($0) * Float(xTB_ZJPerHartree)
    }
  }
  
  func getOrbitalOccupations() {
    let orbitalCount = calculator.orbitals.count
    let orbitalOccupations64 = getDoubleArray(
      xtb_getOrbitalOccupations, size: orbitalCount)
    orbitalOccupations = orbitalOccupations64.map(Float.init)
  }
  
  func getOrbitalCoefficients() {
    let orbitalCount = calculator.orbitals.count
    let orbitalCoefficients64 = getDoubleArray(
      xtb_getOrbitalCoefficients, size: orbitalCount * orbitalCount)
    orbitalCoefficients = orbitalCoefficients64.map(Float.init)
  }
}
