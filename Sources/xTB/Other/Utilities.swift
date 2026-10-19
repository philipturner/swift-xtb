//
//  Utilities.swift
//  swift-xtb
//
//  Created by Philip Turner on 5/30/24.
//

/// Utility function for converting a position array (in nm) to a packed
/// array (in Bohr).
func convertPositions(_ input: [SIMD3<Double>]) -> [Double] {
  var output: [Double] = []
  for atomID in input.indices {
    let positionInNm = input[atomID]
    let positionInBohr = positionInNm * xTB_BohrPerNm
    
    // Add to the packed array.
    for laneID in 0..<3 {
      let element = positionInBohr[laneID]
      output.append(element)
    }
  }
  return output
}

/// Utility function for converting a packed array (in Ha/Bohr) to a force
/// array (in zJ/nm). In addition, converts a gradient into a force by negating
/// it.
func convertGradientToForces(_ input: [Double]) -> [SIMD3<Double>] {
  // The caller should have guaranteed that the input array's length is
  // divisible by 3.
  var output: [SIMD3<Double>] = []
  var cursor: Int = .zero
  for _ in 0..<(input.count / 3) {
    var gradientInHaPerBohr: SIMD3<Double> = .zero
    for laneID in 0..<3 {
      let element = input[cursor]
      gradientInHaPerBohr[laneID] = element
      cursor += 1
    }
    
    let scaleFactor = xTB_ZJPerHartree / xTB_NmPerBohr
    let gradientInZJPerNm = gradientInHaPerBohr * scaleFactor
    let forceInZJPerNm = -gradientInZJPerNm
    output.append(forceInZJPerNm)
  }
  return output
}
