//
//  xTB_Calculator+Caching.swift
//  swift-xtb
//
//  Created by Philip Turner on 6/27/25.
//

import C_xTB

extension xTB_Calculator {
  func invalidateSinglepoint() {
    results = nil
  }
  
  func requestSinglepoint() {
    guard results == nil else {
      return
    }
    
    // Do not access via the public API, otherwise it calls
    // 'invalidateSinglepoint()'. The outcome would be harmless, but we still
    // don't want it.
    //
    // An alternative option would be managing '.molecule' state changes in
    // a separate container, outside of the molecule data types. This removes
    // encapsulation/abstraction and has proven unworkable.
    storage.molecule.update()
    
    let results = xTB_Results()
    xtb_singlepoint(
      xTB_Environment.tEnvironment,
      tMolecule,
      tCalculator,
      results.tResults)
    results.calculator = self
    
    self.results = results
  }
  
  func ensureEnergyCached() {
    requestSinglepoint()
    
    if results.energy == nil {
      results.energy = results.getEnergy()
    }
  }
  
  func ensureMoleculeCached() {
    requestSinglepoint()
    
    if results.forces == nil {
      results.forces = results.getForces()
      
      switch hamiltonian {
      case .forceField:
        results.charges = []
        results.bondOrders = []
      case .tightBinding:
        results.charges = results.getCharges()
        results.bondOrders = results.getBondOrders()
      }
    }
  }
  
  func ensureOrbitalsCached() {
    requestSinglepoint()
    
    if results.orbitalEigenvalues == nil {
      results.checkOrbitalCount()
      
      switch hamiltonian {
      case .forceField:
        results.orbitalEigenvalues = []
        results.orbitalOccupations = []
        results.orbitalCoefficients = []
      case .tightBinding:
        results.orbitalEigenvalues = results.getOrbitalEigenvalues()
        results.orbitalOccupations = results.getOrbitalOccupations()
        results.orbitalCoefficients = results.getOrbitalCoefficients()
      }
    }
  }
}
