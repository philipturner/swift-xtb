//
//  xTB_ExternalCharges.swift
//
//
//  Created by Philip Turner on 5/30/24.
//

/// External charge potential.
public struct xTB_ExternalCharge {
  /// Partial charge in units of proton charge.
  public var charge: Float = .zero
  
  /// Chemical hardness (in atomic units).
  public var chemicalHardness: UInt32 = .zero
  
  /// Coordinates (in nm).
  public var position: SIMD3<Float> = .zero
  
  public init() {
    
  }
}

extension xTB_Calculator {
  /// External charge potential.
  public var externalCharges: [xTB_ExternalCharge] {
    _read {
      yield state.externalCharges
    }
    _modify {
      yield &state.externalCharges
      updateRecord.externalCharges = true
    }
  }
  
  func setExternalCharges(_ externalCharges: [xTB_ExternalCharge]) {
    // Erase the previous external potential.
    xtb_releaseExternalCharges(
      environment.pointer, self.pointer)
    
    // Determine the positions.
    var positions64: [Double] = []
    for chargeID in externalCharges.indices {
      // Convert the position from nm to Bohr.
      let positionInNm = externalCharges[chargeID].position
      let positionInBohr = positionInNm * Float(xTB_BohrPerNm)
      
      // Add to the packed array.
      for laneID in 0..<3 {
        let element = positionInBohr[laneID]
        positions64.append(Double(element))
      }
    }
    
    // Determine the remaining parameters.
    var n = Int32(externalCharges.count)
    var numbers: [Int32] = []
    var charges: [Double] = []
    for chargeID in externalCharges.indices {
      let externalCharge = externalCharges[chargeID]
      numbers.append(Int32(externalCharge.chemicalHardness))
      charges.append(Double(externalCharge.charge))
    }
    
    // Initialize the external potential.
    xtb_setExternalCharges(
      environment.pointer,
      self.pointer,
      &n,
      &numbers,
      &charges,
      &positions64)
  }
}
