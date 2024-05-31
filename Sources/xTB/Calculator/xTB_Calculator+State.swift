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
    var externalCharges: xTB_ExternalCharges!
    var molecule: xTB_Molecule!
    var orbitals: xTB_Orbitals!
  }
  
  struct UpdateRecord {
    var externalCharges: Bool = false
    var molecule: Bool = false
    
    mutating func erase() {
      externalCharges = false
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
    if updateRecord.externalCharges {
      externalCharges.update()
    }
    if updateRecord.molecule {
      molecule.update()
    }
    updateRecord.erase()
  }
  
  /// Run a self-consistent field calculation.
  private func singlepoint() {
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
  func ensureMoleculeCached() {
    
  }
  
  func ensureOrbitalsCached() {
    requestSinglepoint()
    
    if results.orbitalEigenvalues == nil {
      results.getOrbitalEigenvalues()
      results.getOrbitalOccupations()
      results.getOrbitalCoefficients()
    }
  }
}
