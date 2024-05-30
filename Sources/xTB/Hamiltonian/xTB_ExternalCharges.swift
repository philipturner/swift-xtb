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
  public func setExternalCharges(_ externalCharges: [xTB_ExternalCharge]) {
    // Release any previous charges.
    xtb_releaseExternalCharges(environment.pointer, pointer)
    
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
    
    /*
     xtb_setExternalCharges(xtb_TEnvironment /* env */,
                            xtb_TCalculator /* calc */,
                            int* /* n */,
                            int* /* numbers [n] */,
                            double* /* charges [n] */,
                            double* /* positions [n][3] */) XTB_API_SUFFIX__VERSION_1_0_0;
     */
  }
}
