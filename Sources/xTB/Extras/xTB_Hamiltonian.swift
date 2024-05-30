//
//  xTB_Hamiltonian.swift
//
//
//  Created by Philip Turner on 5/30/24.
//

public enum xTB_Hamiltonian {
  /// GFN-FF
  case forceField
  
  /// GFN(n)-xTB
  case tightBinding(Int)
}

extension xTB_Calculator {
  func loadHamiltonian(_ hamiltonian: xTB_Hamiltonian) {
    switch hamiltonian {
    case .forceField:
      xtb_loadGFNFF(
        environment.env, molecule.mol, calc, nil)
    case .tightBinding(let version):
      if version == 0 {
        xtb_loadGFN0xTB(
          environment.env, molecule.mol, calc, nil)
      } else if version == 1 {
        xtb_loadGFN1xTB(
          environment.env, molecule.mol, calc, nil)
      } else if version == 2 {
        xtb_loadGFN2xTB(
          environment.env, molecule.mol, calc, nil)
      } else {
        fatalError("Unsupported tight binding version.")
      }
    }
  }
}

/// Load GFN0-xTB calculator
let xtb_loadGFN0xTB: @convention(c) (
  xtb_TEnvironment?,
  xtb_TMolecule?,
  xtb_TCalculator?,
  UnsafePointer<CChar>? // filename
) -> Void =
xTB_Library.loadSymbol(name: "xtb_loadGFN0xTB")

/// Load GFN1-xTB calculator
let xtb_loadGFN1xTB: @convention(c) (
  xtb_TEnvironment?,
  xtb_TMolecule?,
  xtb_TCalculator?,
  UnsafePointer<CChar>? // filename
) -> Void =
xTB_Library.loadSymbol(name: "xtb_loadGFN1xTB")

/// Load GFN2-xTB calculator
let xtb_loadGFN2xTB: @convention(c) (
  xtb_TEnvironment?,
  xtb_TMolecule?,
  xtb_TCalculator?,
  UnsafePointer<CChar>? // filename
) -> Void =
xTB_Library.loadSymbol(name: "xtb_loadGFN2xTB")

/// Load GFN-FF calculator
let xtb_loadGFNFF: @convention(c) (
  xtb_TEnvironment?,
  xtb_TMolecule?,
  xtb_TCalculator?,
  UnsafePointer<CChar>? // filename
) -> Void =
xTB_Library.loadSymbol(name: "xtb_loadGFNFF")
