//
//  SelfConsistentField.swift
//
//
//  Created by Philip Turner on 2/19/24.
//

/// A configuration for a self-consistent field.
public struct SelfConsistentFieldDescriptor {
  /// Required. The number of protons in each atom's nucleus.
  public var atomicNumbers: [UInt8]?
  
  /// Required. The net charge on each atom during initialization.
  public var charges: [Float]?
  
  /// Required. The minimum number of fragments to split each electron into.
  public var minimumFragmentCount: Int?
  
  /// Required. The position (in nanometers) of each atom's nucleus.
  public var positions: [SIMD3<Float>]?
  
  // We might be able to support specification of bonds. This would work around
  // the issue where odd-numbered atoms necessitate use of UHF. Unsure how such
  // an API would work though. Alternatively, run UHF and switch to RHF once
  // the initial guess has converged. That doesn't seem workable either, though.
  //
  // Solution: the user must specify every other odd-numbered atom as having +
  // or - charge. We should figure out how to document this "best practice" for
  // usage.
  
  /// Required. The net spin on each atom during initialization.
  ///
  /// This value must be divisible by 1/2. If an atom contains an odd number of
  /// electrons, the value must be nonzero.
  public var spins: [Float]?
}

public struct SelfConsistentField {
  // Treat the electrons with UHF for now.
  var spinDownWaveFunctions: [WaveFunction] = []
  var spinUpWaveFunctions: [WaveFunction] = []
  
  public init(descriptor: SelfConsistentFieldDescriptor) {
    guard let atomicNumbers = descriptor.atomicNumbers,
          let charges = descriptor.charges,
          let positions = descriptor.positions,
          let spins = descriptor.spins else {
      fatalError("Descriptor was invalid.")
    }
    guard atomicNumbers.count == charges.count,
          atomicNumbers.count == positions.count,
          atomicNumbers.count == spins.count else {
      fatalError("Descriptor was invalid.")
    }
    
    var minimumPosition = SIMD3<Float>(repeating: .greatestFiniteMagnitude)
    var maximumPosition = SIMD3<Float>(repeating: -.greatestFiniteMagnitude)
    if atomicNumbers.count == .zero {
      // If the software is engineered well, it shouldn't crash when the
      // problem size is zero.
      minimumPosition = .zero
      maximumPosition = .zero
    } else {
      for position in positions {
        maximumPosition.replace(
          with: position, where: position .> maximumPosition)
        minimumPosition.replace(
          with: position, where: position .< minimumPosition)
      }
    }
    
    // Add 20 Bohr of padding (4x the maximum atomic radius) around each atom.
    // Then, round up the world bounds to a power of 2. This will be the cutoff
    // where the multipole expansion is evaluated.
    let origin = (maximumPosition + minimumPosition) / 2
    var size = (maximumPosition - minimumPosition).max()
    size = max(20 + size + 20, 2 * size)
    if size != size.binade {
      size = 2 * size.binade
    }
    
    var octreeDesc = OctreeDescriptor()
    octreeDesc.origin = origin
    octreeDesc.size = size
    let octree = Octree(descriptor: octreeDesc)
    
    for atomID in atomicNumbers.indices {
      let charge = charges[atomID]
      guard charge == charge.rounded(.down) else {
        fatalError("Charge must be an integer.")
      }
      let Z = Int(atomicNumbers[atomID])
      let occupationZ = Z - Int(charge)
      guard occupationZ >= 0 else {
        fatalError("Nucleus had invalid charge.")
      }
      
      // Giving noble gas atoms a negative charge will cause a crash here.
      let occupations = createShellOccupations(Z: occupationZ)
      let effectiveCharges = createEffectiveCharges(Z: Z)
      guard spins[atomID] == 0 else {
        fatalError("Nonzero spins not supported yet.")
      }
      
      // Iterate over the possible quantum numbers for each shell.
      for shellID in occupations.indices {
        let effectiveCharge = effectiveCharges[shellID]
        let occupation = occupations[shellID]
      }
    }
    
    // Test all of this code once you initialize the orbitals. Then, proceed
    // with the remainder of the Hamiltonian construction.
  }
}
