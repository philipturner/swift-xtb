//
//  Ansatz.swift
//
//
//  Created by Philip Turner on 2/19/24.
//

/// A configuration for an initial guess.
public struct AnsatzDescriptor {
  /// Required. The number of protons in each atom's nucleus.
  public var atomicNumbers: [UInt8]?
  
  /// Required. The minimum number of fragments to split each electron into.
  public var minimumFragmentCount: Int?
  
  /// Required. The net charge on each atom.
  public var netCharges: [Int]?
  
  /// Required. Twice the net spin on each atom.
  public var netSpins: [Int]?
  
  /// Required. An octree that specifies the world bounds and minimum basis set.
  public var octree: Octree?
  
  /// Required. The position (in nanometers) of each atom's nucleus.
  public var positions: [SIMD3<Float>]?
}

/// An initial guess to the electronic structure.
public struct Ansatz {
  // Treat the electrons with ROHF. Orthogonalize all of the wavefunctions
  // against each other, except omit the "orthogonalization force" between
  // spin-down and spin-up electrons. When spin-polarized wavefunctions act on
  // spin-neutral ones, halve the "orthogonalization force".
  //
  // Do something similar when iteratively diagonalizing the subspace. Repeat a
  // loop of matrix multiplications and screened orthogonalizations.
  public var spinDownWaveFunctions: [WaveFunction] = []
  public var spinNeutralWaveFunctions: [WaveFunction] = []
  public var spinUpWaveFunctions: [WaveFunction] = []
  
  public init(descriptor: AnsatzDescriptor) {
    // Test this before proceeding with Hamiltonian construction:
    // - Test the expectation values of atomic radii for each shell in each
    //   group IV atom.
    // - Test nonzero spins on N and Cr atoms.
    guard let atomicNumbers = descriptor.atomicNumbers,
          let minimumFragmentCount = descriptor.minimumFragmentCount,
          let netCharges = descriptor.netCharges,
          let netSpins = descriptor.netSpins,
          let octree = descriptor.octree,
          let positions = descriptor.positions else {
      fatalError("Descriptor was invalid.")
    }
    guard atomicNumbers.count == netCharges.count,
          atomicNumbers.count == positions.count,
          atomicNumbers.count == netSpins.count else {
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
    
    #if false
    for atomID in atomicNumbers.indices {
      let charge = charges[atomID]
      let Z = Int(atomicNumbers[atomID])
      let occupationZ = Z - charge
      guard occupationZ >= 0 else {
        fatalError("Nucleus had invalid charge.")
      }
      
      // Giving noble gas atoms a negative charge will cause a crash here.
      // TODO: Change this so it won't crash.
      let occupations = createShellOccupations(Z: occupationZ)
      let effectiveCharges = createEffectiveCharges(Z: Z)
      guard spins[atomID] == 0 else {
        fatalError("Nonzero spins not supported yet.")
      }
      
      // Iterate over the possible quantum numbers for each shell.
      // TODO: Coalesce overlapping spin channels within the electron shell.
      // If there's electrons left over, loop over each spin channel separately.
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
    
    #endif
  }
}
