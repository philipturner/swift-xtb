//
//  SolverTests+Laplacian.swift
//
//
//  Created by Philip Turner on 5/14/24.
//

// Use Neumann boundary conditions for these tests. Also, use the point
// charge model. This would be the multipole expansion of the charge
// distribution created by spreading the nucleus across 8 cells.
//
// NOTE: The Laplacian times the potential does not generate the charge
// density. Replace the right-hand side of the equation with -4πρ.
extension SolverTests {
  static let gridSize: Int = 8
  static let h: Float = 0.25
  
  // The problem size is the number of cells, plus 6 variables for boundary
  // conditions imposed on each cell. To align the matrix rows to the CPU
  // vector width, we pad the number 6 to 8.
  static var cellCount: Int { gridSize * gridSize * gridSize }
  static var n: Int { cellCount + 8 }
  
  enum BoundaryType {
    case dirichlet
    case neumann
  }
  
  // Dirichlet:
  // - Returns the potential evaluated at the center of the boundary face.
  //
  // Neumann:
  // - Set up the Neumann boundaries, normalize to obey Gauss's Law.
  // - Returns an array of fluxes that must be present at the boundary.
  static func createBoundaryConditions(
    type boundaryType: BoundaryType
  ) -> [SIMD8<Float>] {
    // Create an array that represents the boundary values in each cell.
    //
    // Elements of the condition data structure:
    // - [0] = lower X face
    // - [1] = upper X face
    // - [2] = lower Y face
    // - [3] = upper Y face
    // - [4] = lower Z face
    // - [5] = upper Z face
    var conditionGrid = [SIMD8<Float>](
      repeating: .zero, count: gridSize * gridSize * gridSize)
    
    // Iterate over all the boundary cells in the grid. Eventually, we will
    // skip some internal cells to save time.
    for indexZ in 0..<gridSize {
      for indexY in 0..<gridSize {
        // Skip some loop iterations to minimize execution time.
        var indicesX: [Int] = []
        if indexY == 0 || indexY == gridSize - 1 ||
            indexZ == 0 || indexZ == gridSize - 1 {
          for indexX in 0..<gridSize {
            indicesX.append(indexX)
          }
        } else {
          indicesX = [0, gridSize - 1]
        }
        
        for indexX in indicesX {
          // Compute the center of the cell.
          let cellCenterX = (Float(indexX) + 0.5) * h
          let cellCenterY = (Float(indexY) + 0.5) * h
          let cellCenterZ = (Float(indexZ) + 0.5) * h
          let cellCenter = SIMD3<Float>(cellCenterX, cellCenterY, cellCenterZ)
          
          // Determine the condition on each face.
          var faceConditions: SIMD8<Float> = .zero
          for faceID in 0..<6 {
            let coordinateID = faceID / 2
            let signID = faceID % 2
            
            // Compute the center of the face.
            var faceCenter = cellCenter
            let coordinateDelta = (signID == 0) ? Float(-0.5) : 0.5
            faceCenter[coordinateID] += coordinateDelta * h
            
            // Place the nucleus at the midpoint of the 3D grid.
            let nucleusPosition = 0.5 * SIMD3(repeating: Float(gridSize) * h)
            
            // Find the distance and direction from the nucleus.
            let rDelta = faceCenter - nucleusPosition
            let distance = (rDelta * rDelta).sum().squareRoot()
            
            if boundaryType == .dirichlet {
              // The potential is always positive.
              let potential = 1 / distance
              faceConditions[faceID] = potential
            } else {
              // The gradient is always negative.
              let gradient = -1 / (distance * distance)
              
              // Create the flux vector.
              let direction = rDelta / distance
              let flux = gradient * direction
              
              // Select one scalar of the flux vector.
              var faceFlux = flux[coordinateID]
              faceFlux *= (signID == 0) ? -1 : 1
              faceConditions[faceID] = faceFlux
            }
          }
          
          // Erase the conditions on interior faces.
          let indices = SIMD3<Int>(indexX, indexY, indexZ)
          for coordinateID in 0..<3 {
            let index = indices[coordinateID]
            if index != 0 {
              faceConditions[coordinateID * 2 + 0] = .zero
            }
            if index != gridSize - 1 {
              faceConditions[coordinateID * 2 + 1] = .zero
            }
          }
          
          // Store the condition data structure to memory.
          var cellID = indexZ * (gridSize * gridSize)
          cellID += indexY * gridSize + indexX
          conditionGrid[cellID] = faceConditions
        }
      }
    }
    
    // Correct to obey Gauss's Law.
    if boundaryType == .neumann {
      // Integrate the fluxes along the domain boundaries.
      var accumulator: Double = .zero
      for cellID in conditionGrid.indices {
        let faceFluxes = conditionGrid[cellID]
        let fluxTerm = faceFluxes.sum()
        let drTerm = h * h
        accumulator += Double(fluxTerm * drTerm)
      }
      let surfaceIntegral = Float(accumulator)
      
      // Rescale to reflect the charge enclosed.
      let chargeEnclosed: Float = 1
      let actual = surfaceIntegral
      let expected = -4 * Float.pi * chargeEnclosed
      let scaleFactor = expected / actual
      
      for cellID in conditionGrid.indices {
        var faceFluxes = conditionGrid[cellID]
        faceFluxes *= scaleFactor
        conditionGrid[cellID] = faceFluxes
      }
    }
    
    // Return the array of flux data structures.
    return conditionGrid
  }
  
  static func createLaplacianMatrix() -> [Float] {
    // Allocate a matrix.
    var matrix = [Float](repeating: .zero, count: n * n)
    
    // Set the eight extraneous variables to the identity. These variables
    // adapt the boundary conditions to the functional form of a
    // matrix operator.
    for constraintID in cellCount..<cellCount + 8 {
      // Fill in a diagonal of the matrix.
      let diagonalAddress = constraintID * n + constraintID
      matrix[diagonalAddress] = 1
    }
    
    // Fetch the boundary conditions.
    let dirichletConditions = createBoundaryConditions(type: .dirichlet)
    let neumannConditions = createBoundaryConditions(type: .neumann)
    
    // Extra equation to determine the shift in V:
    //
    // Σ w_i * (φ_{extrapolated} - φ_{dirichlet}) = V shift
    var constraintLHS = [Float](repeating: .zero, count: cellCount)
    var constraintFluxSum: Double = .zero
    var constraintDirichletSum: Double = .zero
    var constraintWeightSum: Double = .zero
    
    // Fill in the entries of the matrix.
    for indexZ in 0..<gridSize {
      for indexY in 0..<gridSize {
        for indexX in 0..<gridSize {
          let cellIndices = SIMD3<Int>(indexX, indexY, indexZ)
          var cellID = cellIndices[2] * (gridSize * gridSize)
          cellID += cellIndices[1] * gridSize + cellIndices[0]
          
          // Fetch any possible boundary conditions.
          let facePotentials = dirichletConditions[cellID]
          let faceFluxes = neumannConditions[cellID]
          
          // Iterate over the faces.
          var linkedCellCount: Int = .zero
          for faceID in 0..<6 {
            let coordinateID = faceID / 2
            let signID = faceID % 2
            
            // Locate the neighboring cell.
            var neighborIndices = cellIndices
            neighborIndices[coordinateID] += (signID == 0) ? -1 : 1
            if neighborIndices[coordinateID] >= 0,
                neighborIndices[coordinateID] < gridSize {
              var neighborID = neighborIndices.z * (gridSize * gridSize)
              neighborID += neighborIndices.y * gridSize + neighborIndices.x
              
              // Link this variable to the neighbor.
              linkedCellCount += 1
              
              // Assign 1 / h^2 to the linking entry.
              let linkAddress = cellID * n + neighborID
              let linkEntry: Float = 1 / (h * h)
              matrix[linkAddress] = linkEntry
            } else {
              let facePotential = facePotentials[faceID]
              let faceFlux = faceFluxes[faceID]
              
              // Impose a Neumann boundary condition, as there are no cells to
              // fetch data from.
              let neumannConditionAddress = (cellID * n + cellCount) + faceID
              let neumannConditionEntry: Float = faceFlux / h
              matrix[neumannConditionAddress] = neumannConditionEntry
              
              // Check that there's a neighbor on the other side, to include in
              // the quadratic interpolation.
              var neighborIndices = cellIndices
              neighborIndices[coordinateID] += (signID == 0) ? 1 : -1
              guard neighborIndices[coordinateID] >= 0,
                    neighborIndices[coordinateID] < gridSize else {
                fatalError(
                  "Could not extrapolate potential to Dirichlet boundary.")
              }
              var neighborID = neighborIndices.z * (gridSize * gridSize)
              neighborID += neighborIndices.y * gridSize + neighborIndices.x
              
              // Write to the constraint equation.
              let weight = h * h
              constraintLHS[neighborID] += Float(-1.0 / 8) * weight
              constraintLHS[cellID] += Float(9.0 / 8) * weight
              
              
              let fluxTerm = 3.0 / 8 * faceFlux
              let dirichletTerm = -facePotential
              constraintFluxSum += Double(fluxTerm * weight)
              constraintDirichletSum += Double(dirichletTerm * weight)
              constraintWeightSum += Double(weight)
            }
          }
          
          // Write the entry along the diagonal (most often -6 / h^2).
          let diagonalEntry = -Float(linkedCellCount) / (h * h)
          let diagonalAddress = cellID * n + cellID
          matrix[diagonalAddress] = diagonalEntry
        }
      }
    }
    
    // Collect the FP32 and FP64 values into one array.
    var constraintEquation = [Float](repeating: .zero, count: n)
    
    // Normalize the coefficients with a weight attached to them.
    do {
      let normalizationFactor = 1 / Float(constraintWeightSum)
      for cellID in 0..<cellCount {
        var cellTerm = constraintLHS[cellID]
        cellTerm *= normalizationFactor
        constraintEquation[cellID] = cellTerm
      }
      
      var fluxTerm = Float(constraintFluxSum)
      var dirichletTerm = Float(constraintDirichletSum)
      fluxTerm *= normalizationFactor
      dirichletTerm *= normalizationFactor
      constraintEquation[cellCount + 0] = fluxTerm
      constraintEquation[cellCount + 1] = dirichletTerm
    }
    
    // Set the last term to -1, so it can be occupied by the potential shift.
    constraintEquation[cellCount + 7] = -1
    
    // Write the constraint equation to the last row.
    for columnID in 0..<n {
      let rowID = n - 1
      let address = rowID * n + columnID
      matrix[address] = constraintEquation[columnID]
    }
    
    // Render the Laplacian.
    for rowID in 0..<n {
      for columnID in 0..<n {
        let address = rowID * n + columnID
        let entry = matrix[address]
        var repr = String(format: "%.2f", entry)
        
        // Format the text to fit a constant spacing.
        if entry.sign == .plus {
          repr = " " + repr
        }
        if entry.magnitude >= 10 {
          repr.removeLast()
        }
        
        // Render the text.
        func makeGreen<T: StringProtocol>(_ string: T) -> String {
          "\u{1b}[0;32m\(string)\u{1b}[0m"
        }
        func makeYellow<T: StringProtocol>(_ string: T) -> String {
          "\u{1b}[0;33m\(string)\u{1b}[0m"
        }
        func makeCyan<T: StringProtocol>(_ string: T) -> String {
          "\u{1b}[0;36m\(string)\u{1b}[0m"
        }
        if entry != .zero {
          if rowID < Self.cellCount, columnID < cellCount {
            if rowID == columnID {
              repr = makeCyan(repr)
            } else {
              repr = makeGreen(repr)
            }
          } else {
            repr = makeYellow(repr)
          }
        }
        print(repr, terminator: " ")
      }
      print()
    }
    
    // Render the constraint equation.
    print()
    print("constraint equation:")
    for columnID in 0..<n {
      let entry = Float(constraintEquation[columnID])
      var repr = String(format: "%.2f", entry)
      
      // Format the text to fit a constant spacing.
      if entry.sign == .plus {
        repr = " " + repr
      }
      if entry.magnitude >= 10 {
        repr.removeLast()
      }
      print(repr, terminator: " ")
    }
    print()
    
    return matrix
  }
  
  // Might need a separate function because the matrix with the constraint
  // equation(s) might not be diagonally dominant.
  // - You could also just truncate the output, ignoring the constraint
  //   variables.
}

