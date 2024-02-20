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
  public var charges: [Int]?
  
  /// Required. The minimum number of fragments to split each electron into.
  public var minimumFragmentCount: Int?
  
  /// Required. The position (in nanometers) of each atom's nucleus.
  public var positions: [SIMD3<Float>]?
  
  /// Required. Twice the net spin on each atom during initialization.
  ///
  /// This is the number of spin-up electrons minus the number of
  /// spin-down electrons. If an atom contains an odd number of
  /// electrons, the value must be nonzero.
  public var spinMultiplicities: [Int]?
}

public struct SelfConsistentField {
  // Treat the electrons with ROHF. Orthogonalize all of the wavefunctions
  // against each other, except omit the "orthogonalization force" between
  // spin-down and spin-up electrons. Whenever spin-polarized wavefunctions act
  // on spin-neutral ones, halve the "orthogonalization force".
  //
  // Do something similar when iteratively diagonalizing the subspace. Repeat a
  // loop of matrix multiplications and screened orthogonalizations.
  var spinDownWaveFunctions: [WaveFunction] = []
  var spinNeutralWaveFunctions: [WaveFunction] = []
  var spinUpWaveFunctions: [WaveFunction] = []
  
  public init(descriptor: SelfConsistentFieldDescriptor) {
    guard let atomicNumbers = descriptor.atomicNumbers,
          let charges = descriptor.charges,
          let minimumFragmentCount = descriptor.minimumFragmentCount,
          let positions = descriptor.positions,
          let spinMultiplicities = descriptor.spinMultiplicities else {
      fatalError("Descriptor was invalid.")
    }
    guard atomicNumbers.count == charges.count,
          atomicNumbers.count == positions.count,
          atomicNumbers.count == spinMultiplicities.count else {
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
    
    #if false
    for atomID in atomicNumbers.indices {
      let charge = charges[atomID]
      let Z = Int(atomicNumbers[atomID])
      let occupationZ = Z - charge
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
      for n in occupations.indices {
        let effectiveCharge = effectiveCharges[n]
        var occupation = occupations[n]
        var l = 0
        while l < n && occupation > 0 {
          var m = -l
          while m <= l && occupation > 0 {
            var orbitalDesc = AtomicOrbitalDescriptor()
            orbitalDesc.Z = effectiveCharge
            orbitalDesc.n = n
            orbitalDesc.l = l
            orbitalDesc.m = m
            let atomicOrbital = AtomicOrbital(descriptor: orbitalDesc)
            
            var waveFunctionDesc = WaveFunctionDescriptor()
            waveFunctionDesc.atomicOrbital = atomicOrbital
            waveFunctionDesc.minimumFragmentCount = minimumFragmentCount
            waveFunctionDesc.nucleusPosition = positions[atomID]
            waveFunctionDesc.octree = octree
            let waveFunction = WaveFunction(descriptor: waveFunctionDesc)
            
            if occupation > 0 {
              spinDownWaveFunctions.append(waveFunction)
              occupation -= 1
            }
            if occupation > 0 {
              spinUpWaveFunctions.append(waveFunction)
              occupation -= 1
            }
            
            m += 1
          }
          l += 1
        }
        if occupation != 0 {
          fatalError("Electron shell had unexpected occupation.")
        }
      }
    }
    
    // Check that the system's charge and spin match the sum of the atomic
    // values.
    var netCharge: Int = 0
    var netSpinMultiplicity: Float = 0
    for _ in spinDownWaveFunctions.indices {
      netCharge -= 1
      netSpin -= 1 / 2
    }
    for _ in spinUpWaveFunctions.indices {
      netCharge += 1
      netSpin += 1 / 2
    }
    for charge in charges {
      netCharge += charge
      ne
    }
    
    // Test all of this code once you initialize the orbitals. Then:
    // - Modify createShellOccupations to operate on one spin channel.
    // - Make nonzero spins valid, test on N and Cr atoms.
    // - Proceed with the remainder of the Hamiltonian construction.
    #endif
  }
}
