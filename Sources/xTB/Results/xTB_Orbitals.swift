//
//  xTB_Orbitals.swift
//
//
//  Created by Philip Turner on 5/30/24.
//

public struct xTB_Orbitals {
  unowned var calculator: xTB_Calculator?
  
  // TODO: Check that the retrieved AO count matches this analytical value,
  // otherwise crash.
  public let count: Int
  
  init(descriptor: xTB_CalculatorDescriptor) {
    count = 0
  }
}
