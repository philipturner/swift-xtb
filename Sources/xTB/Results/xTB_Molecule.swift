//
//  xTB_Molecule.swift
//
//
//  Created by Philip Turner on 5/29/24.
//

/// Molecular structure data class.
public struct xTB_Molecule {
  unowned var calculator: xTB_Calculator!
  
  public let atomicNumbers: [UInt8]
  var netCharge: Float
  var netSpin: Float
  var _positions: [SIMD3<Float>] = []
  
  /// Create new molecular structure data
  init(descriptor: xTB_CalculatorDescriptor) {
    guard let atomicNumbers = descriptor.atomicNumbers else {
      fatalError("Descriptor was incomplete.")
    }
    self.atomicNumbers = atomicNumbers
    self.netCharge = descriptor.netCharge
    self.netSpin = descriptor.netSpin
    self._positions = descriptor.positions ?? xTB_Molecule
      .createInitialPositions(atomCount: atomicNumbers.count)
  }
  
  /// Initialization procedure that bypasses an error with atoms having
  /// equal positions.
  static func createInitialPositions(atomCount: Int) -> [SIMD3<Float>] {
    var output: [SIMD3<Float>] = []
    for atomID in 0..<atomCount {
      let scalar = Float(atomID) * 1e-5
      let position = SIMD3(repeating: scalar)
      output.append(position)
    }
    return output
  }
  
  /// Create the reference-counted object from the C API.
  static func createObject(_ molecule: xTB_Molecule) -> xtb_TMolecule {
    // Convert the positions.
    let atomicNumbers = molecule.atomicNumbers
    let positions = molecule._positions
    guard positions.count == atomicNumbers.count else {
      fatalError("Position count did not match atom count.")
    }
    var positions64 = convertPositions(positions)
    
    // Determine the unpaired electron count.
    let unpairedElectronCount = molecule.netSpin.magnitude * 2
    guard var uhf = Int32(exactly: unpairedElectronCount) else {
      fatalError("Net spin must be divisible by 0.5.")
    }
    
    // Create the object.
    var natoms = Int32(atomicNumbers.count)
    var numbers = atomicNumbers.map(Int32.init)
    var charge = Double(molecule.netCharge)
    let mol = xtb_newMolecule(
      xTB_Environment._environment,
      &natoms,
      &numbers,
      &positions64,
      &charge,
      &uhf,
      nil,
      nil)
    guard let mol else {
      fatalError("Could not create new xTB_Molecule.")
    }
    return mol
  }
}

extension xTB_Molecule {
  /// The position of each atom's nucleus (in nanometers).
  public var positions: [SIMD3<Float>] {
    _read {
      yield _positions
    }
    _modify {
      yield &_positions
    }
  }
  
  func update() {
    guard _positions.count == atomicNumbers.count else {
      fatalError("Position count did not match atom count.")
    }
    let positions64 = convertPositions(_positions)
    
    // Update the molecular structure data.
    xtb_updateMolecule(
      xTB_Environment._environment,
      calculator._molecule,
      positions64,
      nil)
  }
  
  /// The force on each atom (in piconewtons).
  public var forces: [SIMD3<Float>] {
    calculator.ensureMoleculeCached()
    return calculator.results.forces!
  }
  
  /// Partial charge in units of proton charge.
  public var charges: [Float] {
    calculator.ensureMoleculeCached()
    return calculator.results.charges!
  }
  
  /// Matrix of Wiberg bond orders.
  ///
  /// Dimensions: (atom count) x (atom count)
  public var bondOrders: [Float] {
    calculator.ensureMoleculeCached()
    return calculator.results.bondOrders!
  }
}
