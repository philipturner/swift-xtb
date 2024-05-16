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
  public var fragmentCount: Int?
  
  /// Required. The net charge on each atom.
  public var netCharges: [Int]?
  
  /// Required. Twice the net spin on each atom.
  public var netSpinPolarizations: [Int]?
  
  /// Required. The position (in nanometers) of each atom's nucleus.
  public var positions: [SIMD3<Float>]?
  
  /// Required. The power-2 size of each wavefunction's octree.
  public var sizeExponent: Int?
  
  public init() {
    
  }
}

/// An initial guess to the electronic structure.
public struct Ansatz {
  public var spinDownWaveFunctions: [WaveFunction] = []
  public var spinNeutralWaveFunctions: [WaveFunction] = []
  public var spinUpWaveFunctions: [WaveFunction] = []
  
  public init(descriptor: AnsatzDescriptor) {
    guard let atomicNumbers = descriptor.atomicNumbers,
          let fragmentCount = descriptor.fragmentCount,
          let netCharges = descriptor.netCharges,
          let netSpinPolarizations = descriptor.netSpinPolarizations,
          let positions = descriptor.positions,
          let sizeExponent = descriptor.sizeExponent else {
      fatalError("Descriptor was invalid.")
    }
    guard atomicNumbers.count == netCharges.count,
          atomicNumbers.count == netSpinPolarizations.count,
          atomicNumbers.count == positions.count else {
      fatalError("Descriptor was invalid.")
    }
    
    for atomID in atomicNumbers.indices {
      // Check that charge and spin are valid.
      let Z = atomicNumbers[atomID]
      let charge = netCharges[atomID]
      let spin = netSpinPolarizations[atomID]
      if charge > Z {
        fatalError("Invalid atom charge.")
      }
      
      let electronCount = Int(Z) - charge
      guard electronCount % 2 == spin.magnitude % 2 else {
        fatalError("Spin multiplicity conflicted with electron count.")
      }
      let spinNeutralOrbitals = (electronCount - Int(spin.magnitude)) / 2
      let spinPolarizedOrbitals = Int(spin.magnitude)
      
      var electronCountDown: Int
      var electronCountUp: Int
      if spin > 0 {
        electronCountDown = spinNeutralOrbitals
        electronCountUp = spinNeutralOrbitals + spinPolarizedOrbitals
      } else {
        electronCountDown = spinNeutralOrbitals + spinPolarizedOrbitals
        electronCountUp = spinNeutralOrbitals
      }
      let spinDownOccupations = createShellOccupations(
        electronCount: electronCountDown)
      let spinNeutralOccupations = createShellOccupations(
        electronCount: spinNeutralOrbitals)
      let spinUpOccupations = createShellOccupations(
        electronCount: electronCountUp)
      let shellCharges = createEffectiveCharges(
        Z: Z, spinDownOccupations: spinDownOccupations,
        spinUpOccupations: spinUpOccupations)
      
      // This loop would be a great place to parallelize across multiple CPU
      // cores. However, we are not permitted to use multicore CPU or hardware
      // acceleration yet. Some very important performance investigations may
      // require a context that is entirely single-core CPU.
      for n in shellCharges.indices {
        let spinDownOccupation = spinDownOccupations[n]
        let spinNeutralOccupation = spinNeutralOccupations[n]
        let spinUpOccupation = spinUpOccupations[n]
        var remainingElectrons = max(spinDownOccupation, spinUpOccupation)
        
        var waveFunctions: [WaveFunction] = []
        var l = 0
        while l < n && remainingElectrons > 0 {
          var m = -l
          while m <= l && remainingElectrons > 0 {
            var orbitalDesc = AtomicOrbitalDescriptor()
            orbitalDesc.Z = shellCharges[n]
            orbitalDesc.n = n
            orbitalDesc.l = l
            orbitalDesc.m = m
            let atomicOrbital = AtomicOrbital(descriptor: orbitalDesc)
            
            var waveFunctionDesc = WaveFunctionDescriptor()
            waveFunctionDesc.atomicOrbital = atomicOrbital
            waveFunctionDesc.fragmentCount = fragmentCount
            waveFunctionDesc.nucleusPosition = positions[atomID]
            waveFunctionDesc.sizeExponent = sizeExponent
            let waveFunction = WaveFunction(descriptor: waveFunctionDesc)
            waveFunctions.append(waveFunction)
            
            remainingElectrons -= 1
            m += 1
          }
          l += 1
        }
        guard remainingElectrons == 0 else {
          fatalError("Unexpected number of remaining electrons.")
        }
        spinDownWaveFunctions += Array(
          waveFunctions[spinNeutralOccupation..<spinDownOccupation])
        spinNeutralWaveFunctions += Array(
          waveFunctions[0..<spinNeutralOccupation])
        spinUpWaveFunctions += Array(
          waveFunctions[spinNeutralOccupation..<spinUpOccupation])
      }
    }
  }
}
