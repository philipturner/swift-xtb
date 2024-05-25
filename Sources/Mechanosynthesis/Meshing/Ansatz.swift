//
//  Ansatz.swift
//
//
//  Created by Philip Turner on 2/19/24.
//

import Dispatch

/// A configuration for an initial guess.
public struct AnsatzDescriptor {
  /// Required. The number of protons in the atom's nucleus.
  public var atomicNumber: UInt8?
  
  /// Required. The minimum number of fragments to split the electron into.
  public var fragmentCount: Int = 1000
  
  /// Required. The net charge on the atom.
  public var netCharge: Int = .zero
  
  /// Required. The net spin on the atom.
  ///
  /// This number can be a fraction. For example, the neutral hydrogen atom
  /// would have a spin of 1/2 or -1/2.
  public var netSpin: Float = .zero
  
  /// Required. The position (in nanometers) of the atom's nucleus.
  public var position: SIMD3<Float>?
  
  /// Required. The power-2 size of the wavefunction's octree.
  ///
  /// Note each orbital's octree is centered at the origin, even if the nucleus
  /// is located far from the origin.
  public var sizeExponent: Int?
  
  public init() {
    
  }
}

/// An initial guess at an atom's electronic structure.
public struct Ansatz {
  public var spinDownOrbitals: [HydrogenicOrbital] = []
  public var spinNeutralOrbitals: [HydrogenicOrbital] = []
  public var spinUpOrbitals: [HydrogenicOrbital] = []
  
  public init(descriptor: AnsatzDescriptor) {
    guard let atomicNumber = descriptor.atomicNumber,
          let position = descriptor.position,
          let sizeExponent = descriptor.sizeExponent else {
      fatalError("Descriptor was incomplete.")
    }
    
    // Check that charge and spin are valid.
    guard let spinPolarization = Int(exactly: descriptor.netSpin * 2) else {
      fatalError("Net spin was not divisible by 1/2.")
    }
    guard descriptor.netCharge <= atomicNumber else {
      fatalError("Net charge was greater than nuclear charge.")
    }
    
    let electronCount = Int(atomicNumber) - descriptor.netCharge
    guard electronCount % 2 == spinPolarization.magnitude % 2 else {
      fatalError("Spin multiplicity conflicted with electron count.")
    }
    let polarizedOrbitalCount = Int(spinPolarization.magnitude)
    let neutralOrbitalCount = (electronCount - polarizedOrbitalCount) / 2
    
    var electronCountDown: Int
    var electronCountUp: Int
    if spinPolarization > 0 {
      electronCountDown = neutralOrbitalCount
      electronCountUp = neutralOrbitalCount + polarizedOrbitalCount
    } else {
      electronCountDown = neutralOrbitalCount + polarizedOrbitalCount
      electronCountUp = neutralOrbitalCount
    }
    let spinDownOccupations = createShellOccupations(
      electronCount: electronCountDown)
    let spinNeutralOccupations = createShellOccupations(
      electronCount: neutralOrbitalCount)
    let spinUpOccupations = createShellOccupations(
      electronCount: electronCountUp)
    let shellCharges = createEffectiveCharges(
      Z: atomicNumber,
      spinDownOccupations: spinDownOccupations,
      spinUpOccupations: spinUpOccupations)
    
    // This loop would be a great place to parallelize across multiple CPU
    // cores. However, we are not permitted to use multicore CPU or hardware
    // acceleration yet. Some very important performance investigations may
    // require a context that is entirely single-core CPU.
    //
    // Even the smallest examples (625 fragments, Z=1) take at least 300 μs.
    
    for n in shellCharges.indices {
      let spinDownOccupation = spinDownOccupations[n]
      let spinNeutralOccupation = spinNeutralOccupations[n]
      let spinUpOccupation = spinUpOccupations[n]
      var remainingElectrons = max(spinDownOccupation, spinUpOccupation)
      
      var orbitals: [HydrogenicOrbital] = []
      var l = 0
      while l < n && remainingElectrons > 0 {
        var m = -l
        while m <= l && remainingElectrons > 0 {
          var basisFunctionDesc = BasisFunctionDescriptor()
          basisFunctionDesc.Z = shellCharges[n]
          basisFunctionDesc.n = n
          basisFunctionDesc.l = l
          basisFunctionDesc.m = m
          let basisFunction = BasisFunction(descriptor: basisFunctionDesc)
          
          var orbitalDesc = HydrogenicOrbitalDescriptor()
          orbitalDesc.basisFunction = basisFunction
          orbitalDesc.fragmentCount = descriptor.fragmentCount
          orbitalDesc.nucleusPosition = position
          orbitalDesc.sizeExponent = sizeExponent
          let orbital = HydrogenicOrbital(descriptor: orbitalDesc)
          orbitals.append(orbital)
          
          remainingElectrons -= 1
          m += 1
        }
        l += 1
      }
      guard remainingElectrons == 0 else {
        fatalError("Unexpected number of remaining electrons.")
      }
      spinNeutralOrbitals += Array(orbitals[0..<spinNeutralOccupation])
      spinDownOrbitals += Array(
        orbitals[spinNeutralOccupation..<spinDownOccupation])
      spinUpOrbitals += Array(
        orbitals[spinNeutralOccupation..<spinUpOccupation])
    }
  }
}
