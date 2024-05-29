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
  /// The occupied spin-orbitals.
  public var orbitals: [HydrogenicOrbital] = []
  
  /// The number of electrons in each orbital.
  public var occupations: [Int] = []
  
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
    
    // Find the number of spin-neutral orbitals.
    let electronCount = Int(atomicNumber) - descriptor.netCharge
    guard electronCount % 2 == spinPolarization.magnitude % 2 else {
      fatalError("Spin multiplicity conflicted with electron count.")
    }
    let polarizedOrbitalCount = Int(spinPolarization.magnitude)
    let neutralOrbitalCount = (electronCount - polarizedOrbitalCount) / 2
    let spinNeutralOccupations = createShellOccupations(
      electronCount: neutralOrbitalCount)
    
    // Find the number of spin-polarized orbitals.
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
    let spinUpOccupations = createShellOccupations(
      electronCount: electronCountUp)
    
    // Find the effective charge in each electron shell.
    let shellCharges = createEffectiveCharges(
      Z: atomicNumber,
      spinDownOccupations: spinDownOccupations,
      spinUpOccupations: spinUpOccupations)
    
    // Find the parameters to initialize every orbital.
    var orbitalDescriptors: [HydrogenicOrbitalDescriptor] = []
    for n in shellCharges.indices {
      let neutralOccupation = spinNeutralOccupations[n]
      let downOccupation = spinDownOccupations[n]
      let upOccupation = spinUpOccupations[n]
      
      var remainingOrbitals = max(downOccupation, upOccupation)
      var remainingNeutralOrbitals = neutralOccupation
      var l = 0
      while l < n && remainingOrbitals > 0 {
        var m = -l
        while m <= l && remainingOrbitals > 0 {
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
          orbitalDescriptors.append(orbitalDesc)
          
          if remainingNeutralOrbitals > 0 {
            occupations.append(2)
          } else {
            occupations.append(1)
          }
          
          remainingOrbitals -= 1
          remainingNeutralOrbitals -= 1
          m += 1
        }
        l += 1
      }
    }
    
    // Create the orbitals.
    orbitals = Self.createOrbitals(descriptors: orbitalDescriptors)
  }
  
  // Generate the octrees for the set of orbitals, while eliding unnecessary
  // computations and parallelizing the rest.
  //
  // The matching stage is technically O(n^2), but it takes ~1000x less
  // time than octree generation (where n ≤ 118). The matching algorithm was
  // chosen because of its simplicity.
  static func createOrbitals(
    descriptors: [HydrogenicOrbitalDescriptor]
  ) -> [HydrogenicOrbital] {
    // De-duplicate the descriptors with the same magnetic quantum number.
    let uniqueDescriptors = descriptors.filter {
      $0.basisFunction!.m == -$0.basisFunction!.l
    }
    
    // Initialize the orbitals in parallel.
    var uniqueOrbitals: [HydrogenicOrbital]
    if uniqueDescriptors.count > 1 {
      let serialQueue = DispatchQueue(label: "Ansatz.init(descriptor:)")
      var completedTasks = [HydrogenicOrbital?](
        repeating: nil, count: uniqueDescriptors.count)
      
      // Dispatch each orbital as a separate task.
      DispatchQueue.concurrentPerform(
        iterations: uniqueDescriptors.count
      ) { z in
        let orbitalDesc = uniqueDescriptors[z]
        let orbital = HydrogenicOrbital(descriptor: orbitalDesc)
        serialQueue.sync {
          completedTasks[z] = orbital
        }
      }
      uniqueOrbitals = serialQueue.sync {
        completedTasks.map { $0! }
      }
    } else {
      uniqueOrbitals = uniqueDescriptors.map {
        HydrogenicOrbital(descriptor: $0)
      }
    }
    
    // Iterate over the duplicated descriptors.
    var output: [HydrogenicOrbital] = []
    for orbitalID in descriptors.indices {
      let duplicatedDescriptor = descriptors[orbitalID]
      let duplicatedFunction = duplicatedDescriptor.basisFunction!
      
      // Iterate over the unique descriptors.
      var matchedDescriptorID: Int?
      for uniqueDescriptorID in uniqueDescriptors.indices {
        let uniqueDescriptor = uniqueDescriptors[uniqueDescriptorID]
        let uniqueFunction = uniqueDescriptor.basisFunction!
        
        // Compare the quantum numbers.
        if duplicatedFunction.n == uniqueFunction.n,
           duplicatedFunction.l == uniqueFunction.l {
          matchedDescriptorID = uniqueDescriptorID
          break
        }
      }
      guard let matchedDescriptorID else {
        fatalError("Could not match to a unique descriptor.")
      }
      
      // Overwrite the matched orbital's basis function with your own.
      var copiedOrbital = uniqueOrbitals[matchedDescriptorID]
      copiedOrbital.basisFunction = duplicatedFunction
      output.append(copiedOrbital)
    }
    return output
  }
}
