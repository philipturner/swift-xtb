//
//  xTB_Molecule.swift
//
//
//  Created by Philip Turner on 5/29/24.
//

/// Molecular structure data class.
public struct xTB_Molecule {
  weak var calculator: xTB_Calculator!
  
  public let atomicNumbers: [UInt8]
  public let netCharge: Float
  public let netSpin: Float
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
}

extension xTB_Molecule {
  func update() {
    guard _positions.count == atomicNumbers.count else {
      fatalError("Position count did not match atom count.")
    }
    let positions64 = Self.convertPositions(_positions)
    
    // Update the molecular structure data.
    xtb_updateMolecule(
      xTB_Environment._environment,
      calculator._molecule,
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
  
  /// Create the reference-counted object from the C API.
  static func createObject(_ molecule: xTB_Molecule) -> xtb_TMolecule {
    // Convert the positions.
    let atomicNumbers = molecule.atomicNumbers
    let positions = molecule._positions
    guard positions.count == atomicNumbers.count else {
      fatalError("Position count did not match atom count.")
    }
    var positions64 = Self.convertPositions(positions)
    
    // Determine the unpaired electron count.
    let uhf = Int32(exactly: molecule.netSpin * 2)
    guard var uhf else {
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
