//
//  xTB_Results.swift
//  
//
//  Created by Philip Turner on 5/30/24.
//

class xTB_Results {
  var _results: xtb_TResults
  
  var energy: Double?
  
  var forces: [SIMD3<Float>]?
  
  var externalChargeForces: [SIMD3<Float>]?
  
  var dipole: SIMD3<Float>?
  
  var charges: [Float]?
  
  var bondOrderMatrix: [Float]?
  
  var orbitalEigenvalues: [Float]?
  
  var orbitalOccupations: [Float]?
  
  var orbitalCoefficientMatrix: [Float]?
  
  /// Create new singlepoint results object
  init() {
    guard let res = xtb_newResults() else {
      fatalError("Could not create new xTB_Results.")
    }
    _results = res
  }
  
  /// Delete singlepoint results object
  deinit {
    xtb_delResults(&_results)
  }
}
