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
    
    func active() -> Bool {
      externalCharges || molecule
    }
  }
  
  func flushUpdateRecord() {
    if updateRecord.externalCharges {
      externalCharges.update()
    }
    if updateRecord.molecule {
      molecule.update()
    }
    updateRecord.erase()
  }
  
  func invalidateSinglepoint() {
    results = nil
  }
}
