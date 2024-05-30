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

/// Returns API version as 10000 * major + 100 * minor + 1 * patch
let xtb_getAPIVersion: @convention(c) () -> Int32 =
xTB_Library.loadSymbol(name: "xtb_getAPIVersion")

//===----------------------------------------------------------------------===//
// Calculation environment
//===----------------------------------------------------------------------===//

/// Create new xtb calculation environment object
let xtb_newEnvironment: @convention(c) () -> xtb_TEnvironment? =
xTB_Library.loadSymbol(name: "xtb_newEnvironment")

/// Delete a xtb calculation environment object
let xtb_delEnvironment: @convention(c) (
  UnsafeMutablePointer<xtb_TEnvironment>) -> Void =
xTB_Library.loadSymbol(name: "xtb_delEnvironment")

/// Check current status of calculation environment
let xtb_checkEnvironment: @convention(c) (
  xtb_TEnvironment) -> Int32 =
xTB_Library.loadSymbol(name: "xtb_checkEnvironment")

/// Show and empty error stack
let xtb_showEnvironment: @convention(c) (
  xtb_TEnvironment,
  UnsafePointer<CChar>? // message
) -> Void =
xTB_Library.loadSymbol(name: "xtb_showEnvironment")

/// Set verbosity of calculation output
let xtb_setVerbosity: @convention(c) (
  xtb_TEnvironment,
  Int32 // verbosity
) -> Void =
xTB_Library.loadSymbol(name: "xtb_setVerbosity")

//===----------------------------------------------------------------------===//
// Molecular structure data class
//===----------------------------------------------------------------------===//

/// Create new molecular structure data (quantities in Bohr)
let xtb_newMolecule: @convention(c) (
  xtb_TEnvironment,
  UnsafePointer<Int32>?, // natoms
  UnsafePointer<Int32>?, // numbers [natoms]
  UnsafePointer<Double>?, // positions [natoms][3]
  UnsafePointer<Double>?, // charge in e
  UnsafePointer<Int32>?, // uhf
  UnsafePointer<Double>?, // lattice [3][3]
  UnsafePointer<CBool>? // periodic [3]
) -> xtb_TMolecule? =
xTB_Library.loadSymbol(name: "xtb_newMolecule")

/// Delete molecular structure data
let xtb_delMolecule: @convention(c) (
  UnsafeMutablePointer<xtb_TMolecule>) -> Void =
xTB_Library.loadSymbol(name: "xtb_delMolecule")

/// Update coordinates and lattice parameters (quantities in Bohr)
let xtb_updateMolecule: @convention(c) (
  xtb_TEnvironment,
  xtb_TMolecule,
  UnsafePointer<Double>?, // positions [natoms][3]
  UnsafePointer<Double>? // lattice [3][3]
) -> Void =
xTB_Library.loadSymbol(name: "xtb_updateMolecule")

//===----------------------------------------------------------------------===//
// Singlepoint calculator
//===----------------------------------------------------------------------===//

/// Create new calculator object
let xtb_newCalculator: @convention(c) () -> xtb_TCalculator? =
xTB_Library.loadSymbol(name: "xtb_newCalculator")

/// Delete calculator object
let xtb_delCalculator: @convention(c) (
  UnsafeMutablePointer<xtb_TCalculator>) -> Void =
xTB_Library.loadSymbol(name: "xtb_delCalculator")

/// Add a external charge potential to calculator (only supported in GFN1/2-xTB)
let xtb_setExternalCharges: @convention(c) (
  xtb_TEnvironment,
  xtb_TCalculator,
  UnsafeMutablePointer<Int32>?, // n
  UnsafeMutablePointer<Int32>?, // numbers [n]
  UnsafeMutablePointer<Double>?, // charges [n]
  UnsafeMutablePointer<Double> // positions [n][3]
) -> Void =
xTB_Library.loadSymbol(name: "xtb_setExternalCharges")

/// Unset the external charge potential
let xtb_releaseExternalCharges: @convention(c) (
  xtb_TEnvironment,
  xtb_TCalculator
) -> Void =
xTB_Library.loadSymbol(name: "xtb_releaseExternalCharges")

/// Perform singlepoint calculation
let xtb_singlepoint: @convention(c) (
  xtb_TEnvironment,
  xtb_TMolecule,
  xtb_TCalculator,
  xtb_TResults
) -> Void =
xTB_Library.loadSymbol(name: "xtb_singlepoint")

//===----------------------------------------------------------------------===//
// Methods of `xtb_TResults`.
//===----------------------------------------------------------------------===//

/// Create new singlepoint results object
let xtb_newResults: @convention(c) () -> xtb_TResults? =
xTB_Library.loadSymbol(name: "xtb_newResults")

/// Delete singlepoint results object
let xtb_delResults: @convention(c) (
  UnsafeMutablePointer<xtb_TResults>) -> Void =
xTB_Library.loadSymbol(name: "xtb_delResults")

/// Query singlepoint results object for energy in Hartree
let xtb_getEnergy: @convention(c) (
  xtb_TEnvironment,
  xtb_TResults,
  UnsafeMutablePointer<Double>? // energy
) -> Void =
xTB_Library.loadSymbol(name: "xtb_getEnergy")

/// Query singlepoint results object for gradient in Hartree / Bohr
let xtb_getGradient: @convention(c) (
  xtb_TEnvironment,
  xtb_TResults,
  UnsafeMutablePointer<Double>? // gradient [natoms][3]
) -> Void =
xTB_Library.loadSymbol(name: "xtb_getGradient")

/// Query singlepoint results object for partial charges in e
let xtb_getCharges: @convention(c) (
  xtb_TEnvironment,
  xtb_TResults,
  UnsafeMutablePointer<Double>? // charges [natoms][3]
) -> Void =
xTB_Library.loadSymbol(name: "xtb_getCharges")
