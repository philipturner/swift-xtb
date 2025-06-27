//
//  xTB_Calculator+State.swift
//
//
//  Created by Philip Turner on 5/30/24.
//

extension xTB_Calculator {
  struct State {
    // Immediately synchronized properties.
    var accuracy: Float = 1.0
    var electronicTemperature: Float = 300
    var maximumIterations: Int = 250
    
    // Lazily synchronized properties.
    var molecule: xTB_Molecule!
    var orbitals: xTB_Orbitals!
  }
  
  struct UpdateRecord {
    var molecule: Bool = false
    
    mutating func erase() {
      molecule = false
    }
  }
  
  func invalidateSinglepoint() {
    results = nil
  }
  
  func requestSinglepoint() {
    if results == nil {
      flushUpdateRecord()
      singlepoint()
    }
  }
  
  /// Ensure the C API objects are up to date.
  private func flushUpdateRecord() {
    print("breakpoint - flushUpdateRecord")
    if updateRecord.molecule {
      molecule.update()
    }
    updateRecord.erase()
  }
  
  /// Run a self-consistent field calculation.
  private func singlepoint() {
    print("breakpoint - xtb_singlepoint")
    let results = xTB_Results()
    xtb_singlepoint(
      xTB_Environment._environment,
      _molecule,
      _calculator,
      results._results)
    results.calculator = self
    self.results = results
  }
}

extension xTB_Calculator {
  func ensureEnergyCached() {
    print("breakpoint - ensureEnergyCached")
    requestSinglepoint()
    
    if results.energy == nil {
      results.getEnergy()
    }
  }
  
  func ensureMoleculeCached() {
    print("breakpoint - ensureMoleculeCached")
    requestSinglepoint()
    
    if results.forces == nil {
      results.getForces()
      results.getCharges()
      results.getBondOrders()
    }
  }
  
  func ensureOrbitalsCached() {
    print("breakpoint - ensureOrbitalsCached")
    requestSinglepoint()
    
    if results.orbitalEigenvalues == nil {
      print("breakpoint - results.orbitalEigenvalues == nil")
      results.checkOrbitalCount()
      results.getOrbitalEigenvalues()
      results.getOrbitalOccupations()
      results.getOrbitalCoefficients()
    }
  }
}
