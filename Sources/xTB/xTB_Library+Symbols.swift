//
//  xTB_Library+Symbols.swift
//
//
//  Created by Philip Turner on 5/29/24.
//

// Opaque pointers to Fortran objects.
typealias xtb_TEnvironment = OpaquePointer
typealias xtb_TMolecule = OpaquePointer
typealias xtb_TCalculator = OpaquePointer
typealias xtb_TResults = OpaquePointer

let xtb_getAPIVersion: @convention(c) () -> Int32 =
xTB_Library.loadSymbol(name: "xtb_getAPIVersion")

let xtb_newEnvironment: @convention(c) () -> xtb_TEnvironment? =
xTB_Library.loadSymbol(name: "xtb_newEnvironment")

let xtb_newCalculator: @convention(c) () -> xtb_TCalculator? =
xTB_Library.loadSymbol(name: "xtb_newCalculator")

let xtb_newResults: @convention(c) () -> xtb_TResults? =
xTB_Library.loadSymbol(name: "xtb_newResults")

let xtb_newMolecule: @convention(c) (
  xtb_TEnvironment?,
  UnsafePointer<Int32>?,
  UnsafePointer<Int32>?,
  UnsafePointer<Double>?,
  UnsafePointer<Double>?,
  UnsafePointer<Int32>?,
  UnsafePointer<Double>?,
  UnsafePointer<CBool>?
) -> xtb_TMolecule? =
xTB_Library.loadSymbol(name: "xtb_newMolecule")

let xtb_updateMolecule: @convention(c) (
  xtb_TEnvironment?,
  xtb_TMolecule?,
  UnsafePointer<Double>?,
  UnsafePointer<Double>?
) -> Void =
xTB_Library.loadSymbol(name: "xtb_updateMolecule")

let xtb_checkEnvironment: @convention(c) (
  xtb_TEnvironment?) -> Int32 =
xTB_Library.loadSymbol(name: "xtb_checkEnvironment")

let xtb_showEnvironment: @convention(c) (
  xtb_TEnvironment?,
  UnsafePointer<CChar>?
) -> Void =
xTB_Library.loadSymbol(name: "xtb_showEnvironment")

let xtb_setVerbosity: @convention(c) (
  xtb_TEnvironment?, Int32) -> Void =
xTB_Library.loadSymbol(name: "xtb_setVerbosity")

let xtb_loadGFN2xTB: @convention(c) (
  xtb_TEnvironment?,
  xtb_TMolecule?,
  xtb_TCalculator?,
  UnsafePointer<CChar>?
) -> Void =
xTB_Library.loadSymbol(name: "xtb_loadGFN2xTB")

let xtb_singlepoint: @convention(c) (
  xtb_TEnvironment?,
  xtb_TMolecule?,
  xtb_TCalculator?,
  xtb_TResults?
) -> Void =
xTB_Library.loadSymbol(name: "xtb_singlepoint")

let xtb_getEnergy: @convention(c) (
  xtb_TEnvironment?,
  xtb_TResults?,
  UnsafeMutablePointer<Double>?
) -> Void =
xTB_Library.loadSymbol(name: "xtb_getEnergy")

let xtb_getGradient: @convention(c) (
  xtb_TEnvironment?,
  xtb_TResults?,
  UnsafeMutablePointer<Double>?
) -> Void =
xTB_Library.loadSymbol(name: "xtb_getGradient")

let xtb_getCharges: @convention(c) (
  xtb_TEnvironment?,
  xtb_TResults?,
  UnsafeMutablePointer<Double>?
) -> Void =
xTB_Library.loadSymbol(name: "xtb_getCharges")

let xtb_delResults: @convention(c) (
  UnsafeMutablePointer<xtb_TResults?>?) -> Void =
xTB_Library.loadSymbol(name: "xtb_delResults")

let xtb_delCalculator: @convention(c) (
  UnsafeMutablePointer<xtb_TCalculator?>?) -> Void =
xTB_Library.loadSymbol(name: "xtb_delCalculator")

let xtb_delMolecule: @convention(c) (
  UnsafeMutablePointer<xtb_TMolecule?>?) -> Void =
xTB_Library.loadSymbol(name: "xtb_delMolecule")

let xtb_delEnvironment: @convention(c) (
  UnsafeMutablePointer<xtb_TEnvironment?>?) -> Void =
xTB_Library.loadSymbol(name: "xtb_delEnvironment")
