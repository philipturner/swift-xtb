//
//  xTB_Molecule.swift
//
//
//  Created by Philip Turner on 5/29/24.
//

/// Molecular structure data class.
class xTB_Molecule {
  var pointer: xtb_TMolecule
  
  // Used for allocating force arrays, etc.
  var atomCount: Int
  
  /// Create new molecular structure data
  init(descriptor: xTB_CalculatorDescriptor) {
    guard let atomicNumbers = descriptor.atomicNumbers,
          let environment = descriptor.environment else {
      fatalError("Descriptor was incomplete.")
    }
    self.atomCount = atomicNumbers.count
    
    // Determine the positions.
    var positions64: [Double] = []
    for atomID in 0..<atomCount {
      // Determine this atom's position.
      var positionInBohr: SIMD3<Float>
      if let positions = descriptor.positions {
        // Convert the position from nm to Bohr.
        let positionInNm = positions[atomID]
        positionInBohr = positionInNm * Float(xTB_BohrPerNm)
      } else {
        // Bypass an error with initializing null positions.
        let scalar = Float(atomID) * 1e-5
        positionInBohr = SIMD3(repeating: scalar)
      }
      
      // Add to the packed array.
      for laneID in 0..<3 {
        let element = positionInBohr[laneID]
        positions64.append(Double(element))
      }
    }
    
    // Determine the unpaired electron count.
    let uhf = Int32(exactly: descriptor.netSpin * 2)
    guard var uhf else {
      fatalError("Net spin must be divisible by 0.5.")
    }
    
    // Create the molecule object.
    var natoms = Int32(atomicNumbers.count)
    var numbers = atomicNumbers.map(Int32.init)
    var charge = Double(descriptor.netCharge)
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
}

extension xTB_Calculator {
  /// Update coordinates (in nm).
  public func setPositions(_ positions: [SIMD3<Float>]) {
    guard positions.count == molecule.atomCount else {
      fatalError("Position count must match atom count.")
    }
    
    // Determine the positions.
    var positions64: [Double] = []
    for atomID in 0..<molecule.atomCount {
      // Convert the position from nm to Bohr.
      let positionInNm = positions[atomID]
      let positionInBohr = positionInNm * Float(xTB_BohrPerNm)
      
      // Add to the packed array.
      for laneID in 0..<3 {
        let element = positionInBohr[laneID]
        positions64.append(Double(element))
      }
    }
    
    // Update the molecular structure data.
    xtb_updateMolecule(
      environment.pointer, molecule.pointer, positions64, nil)
  }
}
