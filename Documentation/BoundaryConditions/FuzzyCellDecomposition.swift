//
//  FuzzyCellDecomposition.swift
//
//
//  Created by Philip Turner on 5/16/24.
//

// Create a continuous charge distribution, and summarize it with FCD-MPE.
//
// This test is a work in progress.
#if false
func testFuzzyCellDecomposition() throws {
  // Set up the water molecule.
  // - H-O bond length: 95.84 pm
  // - H-O-H angle: 104.45°
  func createWaterAnsatz() -> Ansatz {
    // Create the oxygen atom.
    var positions: [SIMD3<Float>] = []
    let oxygenPosition: SIMD3<Float> = .zero
    positions.append(oxygenPosition)
    
    // Create the hydrogen atoms.
    for hydrogenID in 0..<2 {
      let angleDegrees = Float(hydrogenID) * 104.45
      let angleRadians = angleDegrees * (Float.pi / 180)
      
      // Create a direction vector from the angle.
      let directionX = Float.cos(angleRadians)
      let directionZ = -Float.sin(angleRadians)
      let direction = SIMD3(directionX, 0, directionZ)
      
      // Determine the hydrogen atom's position.
      let bohrRadiusInPm: Float = 52.9177210903
      let bondLengthInBohr = 95.84 / bohrRadiusInPm
      let hydrogenPosition = oxygenPosition + direction * bondLengthInBohr
      
      // Append the hydrogen to the array.
      positions.append(hydrogenPosition)
    }
    
    // Specify the topology.
    //
    // The nuclear coordinates range from ±1.8 Bohr.
    // The octree spans from ±0.5 * 16.0 Bohr.
    var ansatzDesc = AnsatzDescriptor()
    ansatzDesc.atomicNumbers = [8, 1, 1]
    ansatzDesc.positions = positions
    ansatzDesc.sizeExponent = 4
    
    // Set the quality to 1000 fragments/electron.
    ansatzDesc.fragmentCount = 1000
    
    // Specify the net charges and spins.
    ansatzDesc.netCharges = [0, 0, 0]
    ansatzDesc.netSpinPolarizations = [0, 1, -1]
    return Ansatz(descriptor: ansatzDesc)
  }
  let ansatz = createWaterAnsatz()
  
  // Inspect each of the wavefunctions.
  //
  // What is the bounding box of each hierarchy level?
  // What is the charge enclosed by each hierarchy level?
  func inspect(waveFunction: WaveFunction) {
    let octree = waveFunction.octree
    inspectFragmentCount()
    inspectBoundingBoxes()
    
    // Inspect the bounding boxes.
    func inspectBoundingBoxes() {
      var boundingBoxes: [Float: (SIMD3<Float>, SIMD3<Float>)] = [:]
      for nodeID in octree.nodes.indices {
        // Retrieve the node and its boundaries.
        let node = octree.nodes[nodeID]
        let boundaryMinimum = node.center - node.spacing / 2
        let boundaryMaximum = node.center + node.spacing / 2
        
        // Retrieve the current bounding box for this level.
        var boundingBox: (SIMD3<Float>, SIMD3<Float>)
        if let previousBoundingBox = boundingBoxes[node.spacing] {
          boundingBox = previousBoundingBox
        } else {
          boundingBox = (
            SIMD3<Float>(repeating: .greatestFiniteMagnitude),
            SIMD3<Float>(repeating: -.greatestFiniteMagnitude))
        }
        
        // Update the bounding box.
        boundingBox.0.replace(
          with: boundaryMinimum, where: boundaryMinimum .< boundingBox.0)
        boundingBox.1.replace(
          with: boundaryMaximum, where: boundaryMaximum .> boundingBox.1)
        boundingBoxes[node.spacing] = boundingBox
      }
      
      // Visualize each level.
      for coordinateID in 0..<3 {
        let coordinateNames: [String] = ["x", "y", "z"]
        let coordinateName = coordinateNames[coordinateID]
        let indentation: String = "  -"
        print(indentation, "\(coordinateName):")
        
        for level in boundingBoxes.keys.sorted().reversed() {
          let indentation: String = "    -"
          let actualLevel = level / 4
          print(indentation, terminator: " ")
          print("levels[\(level) -> \(actualLevel)]:", terminator: " ")
          
          let box = boundingBoxes[level]!
          let minimum = box.0[coordinateID]
          let maximum = box.1[coordinateID]
          print(minimum, "<", coordinateName, "<", maximum)
        }
      }
    }
    
    // Inspect the fragment count.
    func inspectFragmentCount() {
      var fragmentCounts: [Float: Int] = [:]
      for nodeID in octree.nodes.indices {
        // Retrieve the node.
        let node = octree.nodes[nodeID]
        
        // Determine the node's fragment count.
        var nodeFragmentCount: Int = .zero
        for branchID in 0..<8
        where node.branchesMask[branchID] == UInt8.max {
          nodeFragmentCount += 8
        }
        
        // Retrieve the current fragment count for this level.
        var fragmentCount: Int
        if let previousFragmentCount = fragmentCounts[node.spacing] {
          fragmentCount = previousFragmentCount
        } else {
          fragmentCount = .zero
        }
        
        // Update the fragment count.
        fragmentCount += nodeFragmentCount
        fragmentCounts[node.spacing] = fragmentCount
      }
      
      // Visualize the total fragment count.
      do {
        let expectedCount = waveFunction.fragmentCount
        let occupiedCount = fragmentCounts.values.reduce(0, +)
        let allocatedCount = waveFunction.cellValues.count * 8
        print("  - total fragment count:")
        print("    - expected:", expectedCount)
        print("    - occupied:", occupiedCount)
        print("    - allocated:", allocatedCount)
      }
      
      // Visualize each level.
      print("  - fragment distribution:")
      for level in fragmentCounts.keys.sorted().reversed() {
        let indentation: String = "    -"
        let actualLevel = level / 4
        print(indentation, terminator: " ")
        print("levels[\(level) -> \(actualLevel)]:", terminator: " ")
        
        let fragmentCount = fragmentCounts[level]!
        print(fragmentCount)
      }
    }
    
    // Inspect the charge distribution.
    func inspectChargeDistribution() {
      var charges: [Float: Float] = [:]
      for nodeID in octree.nodes.indices {
        // Retrieve the node.
        let node = octree.nodes[nodeID]
        
        // Determine the node's charge.
        var nodeCharge: Float = .zero
        for branchID in 0..<8
        where node.branchesMask[branchID] == UInt8.max {
          let Ψ = waveFunction.cellValues[8 * nodeID + branchID]
          let d3r = node.spacing * node.spacing * node.spacing / 64
          let ΨΨ = Ψ * Ψ * d3r
          nodeCharge += ΨΨ.sum()
        }
        
        // Retrieve the current charge for this level.
        var charge: Float
        if let previousCharge = charges[node.spacing] {
          charge = previousCharge
        } else {
          charge = .zero
        }
      }
    }
  }
  
  print()
  print("O:")
  for spinNeutralOrbitalID in 0..<4 {
    print("- orbitals[\(spinNeutralOrbitalID)]:")
    let waveFunction = ansatz.spinNeutralWaveFunctions[spinNeutralOrbitalID]
    inspect(waveFunction: waveFunction)
  }
  
  print()
  print("H(↑):")
  do {
    print("- orbitals[0]:")
    let waveFunction = ansatz.spinUpWaveFunctions[0]
    inspect(waveFunction: waveFunction)
  }
  
  print()
  print("H(↓):")
  do {
    print("- orbitals[0]:")
    let waveFunction = ansatz.spinDownWaveFunctions[0]
    inspect(waveFunction: waveFunction)
  }
}
#endif
