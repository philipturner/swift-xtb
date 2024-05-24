//
//  FuzzyCellDecomposition.swift
//
//
//  Created by Philip Turner on 5/16/24.
//

// Create a continuous charge distribution, and summarize it with FCD-MPE.
//
// This test is a work in progress.
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
    
    // Set the quality to ≥1000 fragments/electron.
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

func testWorkspace() throws {
  typealias Real = Float
  
  func φ(_ r: SIMD3<Real>) -> Real {
    1 / (r * r).sum().squareRoot()
  }
  
  func φ(_ r: SIMD3<Real>, jitter: Real) -> Real {
    var accumulator: Real = .zero
    for deltaX in 0...1 {
      for deltaY in 0...1 {
        for deltaZ in 0...1 {
          let delta = SIMD3(deltaX, deltaY, deltaZ)
          let rPrime = r + (SIMD3<Real>(delta) * 2 - 1) * jitter
          let value = 1 / (rPrime * rPrime).sum().squareRoot()
          accumulator += value
        }
      }
    }
    return accumulator / 8
  }
  
  // Line FD
  // O(h^2) 1, -2, 1
  // O(h^4) -1/12, 4/3, -5/2, 4/3, -1/12
  // O(h^6) 1/90, -3/20, 3/2, -49/18, 3/2, -3/20, 1/90
  // O(h^8) -1/560, 8/315, -1/5, 8/5, -205/72, 8/5, -1/5, 8/315, -1/560
  let coefficients2: [Real] = [
    1.0,
    -2.0,
    1.0
  ]
  let coefficients4: [Real] = [
    -1.0 / 12, 4.0 / 3,
     -5.0 / 2,
     4.0 / 3, -1.0 / 12
  ]
  let coefficients6: [Real] = [
    1.0 / 90, -3.0 / 20, 3.0 / 2,
    -49.0 / 18,
    3.0 / 2, -3.0 / 20, 1.0 / 90
  ]
  let coefficients8: [Real] = [
    -1.0 / 560, 8.0 / 315, -1.0 / 5, 8.0 / 5,
     -205.0 / 72,
     8.0 / 5, -1.0 / 5, 8.0 / 315, -1.0 / 560
  ]
  let coefficientsMehr: [Real] = [
    -24.0 / 6, 2.0 / 6, 1.0 / 6, 0
  ]
  
  // Try this out on arbitrary functions, whose second derivative is not zero.
  let coordinateJump: Real = 2 / 16
  let h: Real = 0.25 / 16
  let jitter: Real = h / 2
  
  // Second order estimate.
  do {
    var accumulator: Real = .zero
    for coefficientID in 0...2 {
      let coordinate = Real(coefficientID - 1)
      let value = φ(
        SIMD3(5 + coordinateJump * coordinate, 7, 5) / 8,
        jitter: jitter)
      print(SIMD3(5 + coordinateJump * coordinate, 7, 5) / 8, value, coefficients2[coefficientID])
      accumulator += coefficients2[coefficientID] * value
    }
    accumulator = accumulator.magnitude
    print("O(h^2):", accumulator / (h * h) / φ(SIMD3(5, 7, 5) / 8))
  }
  
  // Mehrstellen estimate.
  do {
    var accumulator: Real = .zero
    for deltaX in -1...1 {
      for deltaY in -1...1 {
        for deltaZ in -1...1 {
          let delta = SIMD3(deltaX, deltaY, deltaZ)
          let deltaMagnitude = delta
            .replacing(with: delta &* -1, where: delta .< 0)
          let coefficientID = deltaMagnitude.wrappedSum()
          
          let coordinate = SIMD3<Real>(delta)
          let value = φ(
            SIMD3(5, 7, 5) / 8 + SIMD3(coordinateJump * coordinate) / 8,
            jitter: (deltaX == -1) ? jitter: 0)
          accumulator += coefficientsMehr[coefficientID] * value
        }
      }
    }
    accumulator = accumulator.magnitude
    print("Mehr:  ", accumulator)
  }
  
  // Fourth order estimate.
  do {
    var accumulator: Real = .zero
    for coefficientID in 0...4 {
      let coordinate = Real(coefficientID - 2)
      let value = φ(
        SIMD3(5 + coordinateJump * coordinate, 7, 5) / 8,
        jitter: jitter)
      accumulator += coefficients4[coefficientID] * value
    }
    accumulator = accumulator.magnitude
    print("O(h^4):", accumulator / (h * h) / φ(SIMD3(5, 7, 5) / 8))
  }
  
  // Sixth order estimate.
  do {
    var accumulator: Real = .zero
    for coefficientID in 0...6 {
      let coordinate = Real(coefficientID - 3)
      let value = φ(
        SIMD3(5 + coordinateJump * coordinate, 7, 5) / 8,
        jitter: jitter)
      accumulator += coefficients6[coefficientID] * value
    }
    accumulator = accumulator.magnitude
    print("O(h^6):", accumulator / (h * h) / φ(SIMD3(5, 7, 5) / 8))
  }
  
  // Eighth order estimate.
  do {
    var accumulator: Real = .zero
    for coefficientID in 0...8 {
      let coordinate = Real(coefficientID - 4)
      let value = φ(
        SIMD3(5 + coordinateJump * coordinate, 7, 5) / 8,
        jitter: jitter)
      accumulator += coefficients8[coefficientID] * value
    }
    accumulator = accumulator.magnitude
    print("O(h^8):", accumulator / (h * h) / φ(SIMD3(5, 7, 5) / 8))
  }
}
