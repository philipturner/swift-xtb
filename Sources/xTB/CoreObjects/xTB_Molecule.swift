//
//  xTB_Molecule.swift
//
//
//  Created by Philip Turner on 5/29/24.
//

struct xTB_MoleculeDescriptor {
  var atomicNumbers: [UInt8]?
  var environment: xTB_Environment?
  var netCharge: Float?
  var netSpin: Float?
  var positions: [SIMD3<Float>]?
}

/// Molecular structure data class.
public class xTB_Molecule {
  public let atomicNumbers: [UInt8]
  public let netCharge: Float
  public let netSpin: Float
  
  var pointer: xtb_TMolecule
  weak var calculator: xTB_Calculator!
  
  /// Create new molecular structure data
  init(descriptor: xTB_MoleculeDescriptor) {
    guard let atomicNumbers = descriptor.atomicNumbers,
          let environment = descriptor.environment,
          let netCharge = descriptor.netCharge,
          let netSpin = descriptor.netSpin,
          let positions = descriptor.positions else {
      fatalError("Descriptor was incomplete.")
    }
    self.atomicNumbers = atomicNumbers
    self.netCharge = netCharge
    self.netSpin = netSpin
    
    // Determine the positions.
    guard positions.count == atomicNumbers.count else {
      fatalError("Position count did not match atom count.")
    }
    var positions64 = Self.convertPositions(positions)
    
    // Determine the unpaired electron count.
    let uhf = Int32(exactly: netSpin * 2)
    guard var uhf else {
      fatalError("Net spin must be divisible by 0.5.")
    }
    
    // Create the molecule object.
    var natoms = Int32(atomicNumbers.count)
    var numbers = atomicNumbers.map(Int32.init)
    var charge = Double(netCharge)
    let mol = xtb_newMolecule(
      environment.pointer,
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
    self.pointer = mol
  }
  
  /// Delete molecular structure data.
  deinit {
    xtb_delMolecule(&pointer)
  }
  
  /// The position of each atom's nucleus (in nm).
  public var positions: [SIMD3<Float>] {
    _read {
      yield calculator.state.positions
    }
    _modify {
      yield &calculator.state.positions
      calculator.updateRecord.positions = true
    }
  }
}

extension xTB_Molecule {
  func setPositions(_ positions: [SIMD3<Float>]) {
    guard positions.count == atomicNumbers.count else {
      fatalError("Position count did not match atom count.")
    }
    let positions64 = Self.convertPositions(positions)
    
    // Update the molecular structure data.
    xtb_updateMolecule(
      calculator.environment.pointer,
      self.pointer,
      positions64,
      nil)
  }
}

// MARK: - Utilities

extension xTB_Molecule {
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
  
  /// Utility function for casting a Float32 array (in nm) to a Float64
  /// array (in Bohr).
  static func convertPositions(_ input: [SIMD3<Float>]) -> [Double] {
    var output: [Double] = []
    for elementID in input.indices {
      let positionInNm = input[elementID]
      let positionInBohr = positionInNm * Float(xTB_BohrPerNm)
      
      // Add to the packed array.
      for laneID in 0..<3 {
        let element = positionInBohr[laneID]
        output.append(Double(element))
      }
    }
    return output
  }
}
