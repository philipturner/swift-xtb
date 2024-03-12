//
//  BandForm.swift
//
//
//  Created by Philip Turner on 3/11/24.
//

extension Diagonalization {
  // Returns a matrix of reflectors.
  mutating func reduceToBandForm() -> [Float] {
    var currentReflectors = [Float](
      repeating: 0, count: problemSize * problemSize)
    
    // Reduce the matrix to band form, and collect up the reflectors.
    var blockStart: Int = 0
    while blockStart < problemSize - blockSize {
      // Adjust the loop end, to account for the factorization band offset.
      let blockEnd = min(blockStart + blockSize, problemSize - blockSize)
      defer { blockStart += blockSize }
      
      // Load to panel into the cache, isolating mutations to the matrix A.
      var panel = [Float](repeating: 0, count: blockSize * problemSize)
      for rowID in blockStart..<blockEnd {
        for columnID in 0..<problemSize {
          let matrixAddress = rowID * problemSize + columnID
          let panelAddress = (rowID - blockStart) * problemSize + columnID
          panel[panelAddress] = matrix[matrixAddress]
        }
      }
      
      // Allocate cache memory for the reflectors.
      var panelReflectors = [Float](
        repeating: 0, count: blockSize * problemSize)
      
      // Generate the reflectors.
      for reflectorID in blockStart..<blockEnd {
        // Factor starting at an offset from the diagonal.
        let bandOffset = reflectorID + blockSize
        
        // Load the row into the cache.
        var vector = [Float](repeating: 0, count: problemSize)
        for elementID in 0..<problemSize {
          let address = (reflectorID - blockStart) * problemSize + elementID
          vector[elementID] = panel[address]
        }
        let originalVector = vector
        
        // Apply preceding reflectors (from this panel) to the column.
        for previousReflectorID in blockStart..<reflectorID {
          // Load the reflector into the cache.
          var reflector = [Float](repeating: 0, count: problemSize)
          for elementID in 0..<problemSize {
            let address = (
              previousReflectorID - blockStart) * problemSize + elementID
            reflector[elementID] = panelReflectors[address]
          }
          
          // Apply the reflector.
          var dotProduct: Float = .zero
          for elementID in 0..<problemSize {
            dotProduct += reflector[elementID] * vector[elementID]
          }
          for elementID in 0..<problemSize {
            vector[elementID] -= reflector[elementID] * dotProduct
          }
        }
        
        // Zero out the elements above the band offset.
        for elementID in 0..<bandOffset {
          vector[elementID] = 0
        }
        
        // Take the norm of the vector.
        var norm: Float = .zero
        for elementID in 0..<problemSize {
          norm += vector[elementID] * vector[elementID]
        }
        norm.formSquareRoot()
        
        // Predict the normalization factor.
        let oldSubdiagonal = vector[bandOffset]
        let newSubdiagonal = norm * Float((oldSubdiagonal >= 0) ? -1 : 1)
        let tau = (newSubdiagonal - oldSubdiagonal) / newSubdiagonal
        
        // Check for NANs.
        var nanPresent = false
        let epsilon: Float = 2 * .leastNormalMagnitude
        if (newSubdiagonal - oldSubdiagonal).magnitude < epsilon {
          nanPresent = true
        }
        if newSubdiagonal.magnitude < epsilon {
          nanPresent = true
        }
        
        // Modify the vector, turning it into a reflector.
        for elementID in 0..<problemSize {
          var element = vector[elementID]
          if elementID == bandOffset {
            element = 1
          } else {
            element /= oldSubdiagonal - newSubdiagonal
          }
          element *= tau.squareRoot()
          
          if nanPresent {
            element = 0
          }
          vector[elementID] = element
        }
        
        // Store the reflector to the cache.
        for elementID in 0..<problemSize {
          let address = (reflectorID - blockStart) * problemSize + elementID
          panelReflectors[address] = vector[elementID]
        }
      }
      
      // Apply the reflectors to the matrix, from both sides.
      for directionID in 0..<2 {
        for vectorID in 0..<problemSize {
          var vector = [Float](repeating: 0, count: problemSize)
          if directionID == 0 {
            // Load the row into the cache.
            for elementID in 0..<problemSize {
              let address = vectorID * problemSize + elementID
              vector[elementID] = matrix[address]
            }
          } else {
            // Load the column into the cache.
            for elementID in 0..<problemSize {
              let address = elementID * problemSize + vectorID
              vector[elementID] = matrix[address]
            }
          }
          
          for reflectorID in blockStart..<blockEnd {
            // Load the reflector into the cache.
            var reflector = [Float](repeating: 0, count: problemSize)
            for elementID in 0..<problemSize {
              let address = (
                reflectorID - blockStart) * problemSize + elementID
              reflector[elementID] = panelReflectors[address]
            }
            
            // Apply the reflector.
            var dotProduct: Float = .zero
            for elementID in 0..<problemSize {
              dotProduct += reflector[elementID] * vector[elementID]
            }
            for elementID in 0..<problemSize {
              vector[elementID] -= reflector[elementID] * dotProduct
            }
          }
          
          if directionID == 0 {
            // Store the row to main memory.
            for elementID in 0..<problemSize {
              let address = vectorID * problemSize + elementID
              matrix[address] = vector[elementID]
            }
          } else {
            // Store the column to main memory.
            for elementID in 0..<problemSize {
              let address = elementID * problemSize + vectorID
              matrix[address] = vector[elementID]
            }
          }
        }
      }
      
      // Store the reflectors to main memory.
      for reflectorID in blockStart..<blockEnd {
        for elementID in 0..<problemSize {
          let cacheAddress = (
            reflectorID - blockStart) * problemSize + elementID
          let memoryAddress = reflectorID * problemSize + elementID
          currentReflectors[memoryAddress] = panelReflectors[cacheAddress]
        }
      }
    }
    
    return currentReflectors
  }
}
