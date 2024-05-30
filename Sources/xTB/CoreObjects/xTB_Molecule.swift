//
//  xTB_Molecule.swift
//
//
//  Created by Philip Turner on 5/29/24.
//

/// Molecular structure data class.
class xTB_Molecule {
  var mol: xtb_TMolecule?
  
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
    var positions: [SIMD3<Float>]
    if descriptor.positions != nil {
      positions = descriptor.positions!
    } else {
      positions = []
      for atomID in 0..<atomCount {
        // Bypass an error with initializing null positions.
        let scalar = Float(atomID) * 1e-5
        let position = SIMD3<Float>(repeating: scalar)
        positions.append(position)
      }
    }
    var positions64: [Double] = []
    for atomID in 0..<atomCount {
      let position = positions[atomID]
      positions64.append(Double(position[0]))
      positions64.append(Double(position[1]))
      positions64.append(Double(position[2]))
    }
    
    // Create the molecule object.
    var natoms = Int32(atomicNumbers.count)
    var numbers = atomicNumbers.map(Int32.init)
    var charge = Double(descriptor.netCharge)
    guard var uhf = Int32(exactly: descriptor.netSpin) else {
      fatalError("Net spin must be divisible by 0.5.")
    }
    mol = xtb_newMolecule(
      environment.env,
      &natoms,
      &numbers,
      &positions64,
      &charge,
      &uhf,
      nil,
      nil)
    guard mol != nil else {
      fatalError("Could not create new xTB_Molecule.")
    }
  }
  
  /// Delete molecular structure data.
  deinit {
    xtb_delMolecule(&mol)
  }
}

extension xTB_Calculator {
  /// Update coordinates (in nm).
  public func setPositions(_ positions: [SIMD3<Float>]) {
    guard positions.count == molecule.atomCount else {
      fatalError("Position count must match atom count.")
    }
    
    var output = [Double](repeating: .zero, count: positions.count * 3)
    for atomID in positions.indices {
      let positionInNm = positions[atomID]
      let positionInBohr = positionInNm * Float(xTB_BohrPerNm)
      for laneID in 0..<3 {
        let element = positionInBohr[laneID]
        output[atomID &* 3 &+ laneID] = Double(element)
      }
    }
    xtb_updateMolecule(environment.env, molecule.mol, output, nil)
  }
}
