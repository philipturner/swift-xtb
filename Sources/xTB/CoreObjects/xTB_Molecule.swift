//
//  xTB_Molecule.swift
//
//
//  Created by Philip Turner on 5/29/24.
//

/// A configuration for a molecular structure data class.
struct xTB_MoleculeDescriptor {
  var atomicNumbers: [UInt8]?
  
  var environment: xTB_Environment?
  
  var netCharge: Float = .zero
  
  var netSpin: Float = .zero
  
  init() {
    
  }
}

/// Molecular structure data class.
class xTB_Molecule {
  var mol: xtb_TMolecule?
  
  // Keep a reference to the environment, so it is never deallocated while the
  // molecule is still in use.
  var environment: xTB_Environment
  
  // Used for allocating force arrays, etc.
  var atomCount: Int
  
  /// Create new molecular structure data
  init(descriptor: xTB_MoleculeDescriptor) {
    guard let atomicNumbers = descriptor.atomicNumbers,
          let environment = descriptor.environment else {
      fatalError("Descriptor was incomplete.")
    }
    self.environment = environment
    self.atomCount = atomicNumbers.count
    
    var natoms = Int32(atomicNumbers.count)
    var numbers = atomicNumbers.map(Int32.init)
    var positions = [Double](repeating: .zero, count: Int(natoms * 3))
    do {
      // Bypass an initialization error.
      for elementID in positions.indices {
        positions[elementID] = Double(elementID) * 1e-5
      }
    }
    var charge = Double(descriptor.netCharge)
    guard var uhf = Int32(exactly: descriptor.netSpin) else {
      fatalError("Net spin must be divisible by 0.5.")
    }
    
    mol = xtb_newMolecule(
      environment.env,
      &natoms,
      &numbers,
      &positions,
      &charge,
      &uhf,
      nil,
      nil)
    guard mol != nil else {
      print(environment.status)
      environment.show()
      fatalError("Could not create new xTB_Molecule.")
    }
  }
  
  /// Delete molecular structure data.
  deinit {
    xtb_delMolecule(&mol)
  }
  
  /// Update coordinates (in nm).
  func setPositions(_ positions: [SIMD3<Float>]) {
    var output = [Double](repeating: .zero, count: positions.count * 3)
    for atomID in positions.indices {
      let positionInNm = positions[atomID]
      let positionInBohr = positionInNm * Float(xTB_BohrPerNm)
      for laneID in 0..<3 {
        let element = positionInBohr[laneID]
        output[atomID &* 3 &+ laneID] = Double(element)
      }
    }
    xtb_updateMolecule(environment.env, mol, output, nil)
  }
}
