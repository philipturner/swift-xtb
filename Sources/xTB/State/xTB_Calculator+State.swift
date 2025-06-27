//
//  xTB_Calculator+State.swift
//  swift-xtb
//
//  Created by Philip Turner on 5/30/24.
//

import C_xTB

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
  
  func invalidateSinglepoint() {
    results = nil
  }
  
  func requestSinglepoint() {
    if results == nil {
      if positionsUpdated {
        molecule.update()
      }
      positionsUpdated = false
      
      let results = xTB_Results()
      xtb_singlepoint(
        xTB_Environment.tEnvironment,
        tMolecule,
        tCalculator,
        results.tResults)
      results.calculator = self
      
      self.results = results
    }
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
      results.charges = results.getCharges()
      results.getBondOrders()
    }
  }
  
  func ensureOrbitalsCached() {
    requestSinglepoint()
    
    if results.orbitalEigenvalues == nil {
      results.checkOrbitalCount()
      results.getOrbitalEigenvalues()
      results.getOrbitalOccupations()
      results.getOrbitalCoefficients()
    }
  }
}
