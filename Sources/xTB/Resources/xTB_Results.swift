//
//  xTB_Results.swift
//  
//
//  Created by Philip Turner on 5/30/24.
//

class xTB_Results {
  var pointer: xtb_TResults
  
  /// Create new singlepoint results object
  init() {
    guard let res = xtb_newResults() else {
      fatalError("Could not create new xTB_Results.")
    }
    pointer = res
  }
  
  /// Delete singlepoint results object
  deinit {
    xtb_delResults(&pointer)
  }
}
