//
//  xTB_Library+Symbols.swift
//  swift-xtb
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
  UnsafeMutablePointer<xtb_TEnvironment?>?) -> Void =
xTB_Library.loadSymbol(name: "xtb_delEnvironment")

/// Check current status of calculation environment
let xtb_checkEnvironment: @convention(c) (
  xtb_TEnvironment?) -> Int32 =
xTB_Library.loadSymbol(name: "xtb_checkEnvironment")

/// Show and empty error stack
let xtb_showEnvironment: @convention(c) (
  xtb_TEnvironment?,
  UnsafePointer<CChar>? // message
) -> Void =
xTB_Library.loadSymbol(name: "xtb_showEnvironment")

/// Return and empty error stack
let xtb_getError: @convention(c) (
  xtb_TEnvironment?,
  UnsafeMutablePointer<UInt8>?, // buffer
  UnsafePointer<Int32> // buffersize
) -> Void =
xTB_Library.loadSymbol(name: "xtb_getError")

/// Bind output from this environment
let xtb_setOutput: @convention(c) (
  xtb_TEnvironment?,
  UnsafePointer<CChar>? // filename
) -> Void =
xTB_Library.loadSymbol(name: "xtb_setOutput")

/// Set verbosity of calculation output
let xtb_setVerbosity: @convention(c) (
  xtb_TEnvironment?,
  Int32 // verbosity
) -> Void =
xTB_Library.loadSymbol(name: "xtb_setVerbosity")

//===----------------------------------------------------------------------===//
// Molecular structure data class
//===----------------------------------------------------------------------===//

/// Create new molecular structure data (quantities in Bohr)
let xtb_newMolecule: @convention(c) (
  xtb_TEnvironment?,
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
  UnsafeMutablePointer<xtb_TMolecule?>?) -> Void =
xTB_Library.loadSymbol(name: "xtb_delMolecule")

/// Update coordinates and lattice parameters (quantities in Bohr)
let xtb_updateMolecule: @convention(c) (
  xtb_TEnvironment?,
  xtb_TMolecule?,
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
  UnsafeMutablePointer<xtb_TCalculator?>?) -> Void =
xTB_Library.loadSymbol(name: "xtb_delCalculator")

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

/// Set numerical accuracy of calculator in the range of 1000 to 0.0001
let xtb_setAccuracy: @convention(c) (
  xtb_TEnvironment?,
  xtb_TCalculator?,
  Double // accuracy
) -> Void =
xTB_Library.loadSymbol(name: "xtb_setAccuracy")

/// Set maximum number of iterations for self-consistent TB calculators
let xtb_setMaxIter: @convention(c) (
  xtb_TEnvironment?,
  xtb_TCalculator?,
  Int32 // iterations
) -> Void =
xTB_Library.loadSymbol(name: "xtb_setMaxIter")

/// Set electronic temperature for level filling in tight binding calculators in K
let xtb_setElectronicTemp: @convention(c) (
  xtb_TEnvironment?,
  xtb_TCalculator?,
  Double // temperature
) -> Void =
xTB_Library.loadSymbol(name: "xtb_setElectronicTemp")

/// Perform singlepoint calculation
let xtb_singlepoint: @convention(c) (
  xtb_TEnvironment?,
  xtb_TMolecule?,
  xtb_TCalculator?,
  xtb_TResults?
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
  UnsafeMutablePointer<xtb_TResults?>?) -> Void =
xTB_Library.loadSymbol(name: "xtb_delResults")

/// Create copy from a singlepoint results object
let xtb_copyResults: @convention(c) (
  xtb_TResults?) -> xtb_TResults? =
xTB_Library.loadSymbol(name: "xtb_copyResults")

/// Query singlepoint results object for energy in Hartree
let xtb_getEnergy: @convention(c) (
  xtb_TEnvironment?,
  xtb_TResults?,
  UnsafeMutablePointer<Double>? // energy
) -> Void =
xTB_Library.loadSymbol(name: "xtb_getEnergy")

/// Query singlepoint results object for gradient in Hartree / Bohr
let xtb_getGradient: @convention(c) (
  xtb_TEnvironment?,
  xtb_TResults?,
  UnsafeMutablePointer<Double>? // gradient [natoms][3]
) -> Void =
xTB_Library.loadSymbol(name: "xtb_getGradient")

/// Query singlepoint results object for partial charges in e
let xtb_getCharges: @convention(c) (
  xtb_TEnvironment?,
  xtb_TResults?,
  UnsafeMutablePointer<Double>? // charges [natoms][3]
) -> Void =
xTB_Library.loadSymbol(name: "xtb_getCharges")

/// Query singlepoint results object for bond orders
let xtb_getBondOrders: @convention(c) (
  xtb_TEnvironment?,
  xtb_TResults?,
  UnsafeMutablePointer<Double>? // wbo [natoms][natoms]
) -> Void =
xTB_Library.loadSymbol(name: "xtb_getBondOrders")

/// Query singlepoint results object for the number of basis functions
let xtb_getNao: @convention(c) (
  xtb_TEnvironment?,
  xtb_TResults?,
  UnsafeMutablePointer<Int32>? // nao
) -> Void =
xTB_Library.loadSymbol(name: "xtb_getNao")

/// Query singlepoint results object for orbital energies in Hartree [nao]
let xtb_getOrbitalEigenvalues: @convention(c) (
  xtb_TEnvironment?,
  xtb_TResults?,
  UnsafeMutablePointer<Double>? // emo
) -> Void =
xTB_Library.loadSymbol(name: "xtb_getOrbitalEigenvalues")

/// Query singlepoint results object for occupation numbers [nao]
let xtb_getOrbitalOccupations: @convention(c) (
  xtb_TEnvironment?,
  xtb_TResults?,
  UnsafeMutablePointer<Double>? // focc
) -> Void =
xTB_Library.loadSymbol(name: "xtb_getOrbitalOccupations")

/// Query singlepoint results object for orbital coefficients [nao][nao]
let xtb_getOrbitalCoefficients: @convention(c) (
  xtb_TEnvironment?,
  xtb_TResults?,
  UnsafeMutablePointer<Double>? // c
) -> Void =
xTB_Library.loadSymbol(name: "xtb_getOrbitalCoefficients")

//===----------------------------------------------------------------------===//
// SCC eigensolve precision
//
// Only in xTB builds from the 'diagonalization' branch, so these are loaded
// as optionals: nil means the library runs every eigensolve in FP64.
//===----------------------------------------------------------------------===//

/// Set the precision of the SCC eigensolves of an xTB calculator: FP32 while
/// the charges change by more than fp32_rmsq (0: FP64 only), the FP32
/// attempt capped at fp32_maxiter cycles before an FP64 rerun, fp32_refine
/// 1 to refine each FP32 eigensolve in FP64
let xtb_setSCCPrecision: (@convention(c) (
  xtb_TEnvironment?,
  xtb_TCalculator?,
  Double, // fp32_rmsq
  Int32, // fp32_maxiter
  Int32 // fp32_refine
) -> Void)? =
xTB_Library.loadOptionalSymbol(name: "xtb_setSCCPrecision")

/// Query singlepoint results object for the SCC eigensolve precision used and
/// its outcome
let xtb_getSCCPrecision: (@convention(c) (
  xtb_TEnvironment?,
  xtb_TResults?,
  UnsafeMutablePointer<Double>?, // fp32_rmsq
  UnsafeMutablePointer<Int32>?, // fp32_maxiter
  UnsafeMutablePointer<Int32>?, // fp32_refine
  UnsafeMutablePointer<Int32>?, // fp32_cycles
  UnsafeMutablePointer<Int32>?, // fp64_restart
  UnsafeMutablePointer<Int32>? // scc_cycles
) -> Void)? =
xTB_Library.loadOptionalSymbol(name: "xtb_getSCCPrecision")

/// Whether both functions of the SCC precision API are present.
let xtb_supportsSCCPrecision: Bool = {
  switch (xtb_setSCCPrecision != nil, xtb_getSCCPrecision != nil) {
  case (true, true):
    return true
  case (false, false):
    return false
  default:
    fatalError("""
      The xTB library exports only one of xtb_setSCCPrecision and \
      xtb_getSCCPrecision.
      """)
  }
}()
