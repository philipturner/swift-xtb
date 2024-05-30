//
//  xTB_CalculatorStorage.swift
//
//
//  Created by Philip Turner on 5/30/24.
//

/// Stores the state of the calculator.
struct xTB_CalculatorStorage {
  // The backing calculation environment.
  var environment: xTB_Environment
  
  // The backing molecule.
  var molecule: xTB_Molecule
  
  var accuracy: Float = 1.0
  
  var electronicTemperature: Float = 300
  
  var externalCharges: [xTB_ExternalCharge] = []
  
  var maximumIterations: Int = 250
  
  init(environment: xTB_Environment, molecule: xTB_Molecule) {
    self.environment = environment
    self.molecule = molecule
  }
}
