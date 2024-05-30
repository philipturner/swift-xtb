//
//  xTB_Calculator+State.swift
//
//
//  Created by Philip Turner on 5/30/24.
//

extension xTB_Calculator {
  struct State {
    var accuracy: Float = 1.0
    var electronicTemperature: Float = 300
    var externalCharges: [xTB_ExternalCharge] = []
    var maximumIterations: Int = 250
    var positions: [SIMD3<Float>] = []
    
    init() {
      
    }
  }
  
  struct UpdateRecord {
    var accuracy: Bool = false
    var electronicTemperature: Bool = false
    var externalCharges: Bool = false
    var maximumIterations: Bool = false
    var positions: Bool = false
    
    func active() -> Bool {
      accuracy ||
      electronicTemperature ||
      externalCharges ||
      maximumIterations ||
      positions
    }
  }
  
  func flushUpdateRecord() {
    if updateRecord.accuracy {
      setAccuracy(state.accuracy)
    }
    if updateRecord.electronicTemperature {
      setElectronicTemperature(state.electronicTemperature)
    }
    if updateRecord.externalCharges {
      setExternalCharges(state.externalCharges)
    }
    if updateRecord.maximumIterations {
      setMaximumIterations(state.maximumIterations)
    }
    if updateRecord.positions {
      setPositions(state.positions)
    }
    
    // Erase the update record.
    updateRecord = xTB_UpdateRecord()
  }
}
