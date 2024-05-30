//
//  xTB_UpdateRecord.swift
//  
//
//  Created by Philip Turner on 5/30/24.
//

struct xTB_UpdateRecord {
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

extension xTB_Calculator {
  func flushUpdateRecord() {
    if updateRecord.accuracy {
      setAccuracy(accuracy)
    }
    if updateRecord.electronicTemperature {
      setElectronicTemperature(storage.electronicTemperature)
    }
    if updateRecord.externalCharges {
      setExternalCharges(storage.externalCharges)
    }
    if updateRecord.maximumIterations {
      setMaximumIterations(storage.maximumIterations)
    }
    if updateRecord.positions {
      setPositions(storage.molecule.positions)
    }
    
    // Erase the update record.
    updateRecord = xTB_UpdateRecord()
  }
}
