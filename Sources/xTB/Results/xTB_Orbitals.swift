//
//  xTB_Orbitals.swift
//
//
//  Created by Philip Turner on 5/30/24.
//

public struct xTB_Orbitals {
  unowned var calculator: xTB_Calculator!
  
  public let count: Int
  
  init(descriptor: xTB_CalculatorDescriptor) {
    guard let atomicNumbers = descriptor.atomicNumbers else {
      fatalError("Atomic numbers were not specified.")
    }
    if descriptor.hamiltonian == .forceField {
      count = .zero
    } else {
      count = Self.estimateOrbitalCount(atomicNumbers: atomicNumbers)
    }
  }
}

// MARK: - Utilities

extension xTB_Orbitals {
  static func estimateOrbitalCount(atomicNumbers: [UInt8]) -> Int {
    var output: Int = .zero
    for atomicNumber in atomicNumbers {
      guard 0 < atomicNumber, atomicNumber <= 86 else {
        fatalError("Atomic number is out of range.")
      }
      
      // Handle the atoms with a unique basis.
      if atomicNumber == 1 {
        // 1s
        output += 1
        continue
      } else if atomicNumber == 2 {
        // 1s, 2p
        output += 1 + 3
        continue
      } else if atomicNumber == 10 {
        // 2s, 2p, 3d
        output += 1 + 3 + 5
        continue
      }
      
      // Handle the group 1 elements.
      let nobleGases: [UInt8] = [2, 10, 18, 36, 54, 86]
      if nobleGases.contains(atomicNumber - 1) {
        // ns, np
        output += 1 + 3
        continue
      }
      
      // Handle the remaining elements with 4 orbitals.
      if 2 <= atomicNumber, atomicNumber <= 9 {
        // 2s, 2p
        output += 1 + 3
        continue
      } else if atomicNumber == 30 {
        // 4s, 4p
        output += 1 + 3
        continue
      } else if atomicNumber == 48 {
        // 5s, 5p
        output += 1 + 3
        continue
      } else if 80 <= atomicNumber, atomicNumber <= 84 {
        // 6s, 6p
        output += 1 + 3
        continue
      }
      
      // The remaining elements must have 9 orbitals.
      // nspd
      // nd, (n + 1)sp
      output += 1 + 3 + 5
      continue
    }
    return output
  }
}
