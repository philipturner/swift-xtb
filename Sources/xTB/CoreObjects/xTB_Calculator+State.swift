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
    var externalCharges: [xTB_ExternalCharge] = []
    var positions: [SIMD3<Float>] = []
  }
  
  struct UpdateRecord {
    var externalCharges: Bool = false
    var positions: Bool = false
    
    mutating func erase() {
      externalCharges = false
      positions = false
    }
    
    func active() -> Bool {
      externalCharges || positions
    }
  }
  
  func flushUpdateRecord() {
    if updateRecord.externalCharges {
      setExternalCharges(state.externalCharges)
    }
    if updateRecord.positions {
      molecule.setPositions(state.positions)
    }
    updateRecord.erase()
  }
}
