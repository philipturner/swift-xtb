//
//  BandForm.swift
//
//
//  Created by Philip Turner on 3/11/24.
//

import Accelerate

extension Diagonalization {
  // Returns a matrix of reflectors.
  mutating func reduceToBandForm() -> (
    reflectors: [Float], tau: [Float]
  ) {
    // TODO: Accelerate this step with the compact WY transform. Start by
    // generating the T matrix after all reflectors in the panel have been
    // completed.
    var currentReflectors = [Float](
      repeating: 0, count: problemSize * problemSize)
    var currentTau = [Float](repeating: 0, count: problemSize)
    
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
      var panelTau = [Float](repeating: 0, count: blockSize)
      
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
        
        // Apply preceding reflectors (from this panel) to the column.
        for previousReflectorID in blockStart..<reflectorID {
          // Load the reflector into the cache.
          var reflector = [Float](repeating: 0, count: problemSize)
          for elementID in 0..<problemSize {
            let address = (
              previousReflectorID - blockStart) * problemSize + elementID
            reflector[elementID] = panelReflectors[address]
          }
          let tau = panelTau[previousReflectorID - blockStart]
          
          // Apply the reflector.
          var dotProduct: Float = .zero
          for elementID in 0..<problemSize {
            dotProduct += reflector[elementID] * vector[elementID]
          }
          for elementID in 0..<problemSize {
            vector[elementID] -= tau * reflector[elementID] * dotProduct
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
        var tau = (newSubdiagonal - oldSubdiagonal) / newSubdiagonal
        
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
          if elementID < bandOffset {
            element = 0
          } else if elementID == bandOffset {
            element = 1
          } else {
            element /= oldSubdiagonal - newSubdiagonal
          }
          
          if nanPresent {
            element = 0
          }
          vector[elementID] = element
        }
        if nanPresent {
          tau = 0
        }
        
        // Store the reflector to the cache.
        for elementID in 0..<problemSize {
          let address = (reflectorID - blockStart) * problemSize + elementID
          panelReflectors[address] = vector[elementID]
        }
        panelTau[reflectorID - blockStart] = tau
      }
      
      // Perform a GEMM with non-square matrices.
      var reflectorDotProducts = [Float](
        repeating: 0, count: blockSize * blockSize)
      for lhsID in 0..<blockSize {
        for rhsID in 0..<blockSize {
          var dotProduct: Float = .zero
          for elementID in 0..<problemSize {
            let lhsAddress = lhsID * problemSize + elementID
            let rhsAddress = rhsID * problemSize + elementID
            let lhsValue = panelReflectors[lhsAddress]
            let rhsValue = panelReflectors[rhsAddress]
            dotProduct += lhsValue * rhsValue
          }
          let address = lhsID * blockSize + rhsID
          reflectorDotProducts[address] = dotProduct
        }
      }
      
      // Generate the diagonal entries.
      var currentMatrixT = [Float](repeating: 0, count: problemSize * problemSize)
      for transformID in blockStart..<blockEnd {
        let diagonalAddress = transformID * problemSize + transformID
        let dotProductAddress = (
          transformID - blockStart) * blockSize + (transformID - blockStart)
        let dotProduct = reflectorDotProducts[dotProductAddress]
        currentMatrixT[diagonalAddress] = panelTau[transformID - blockStart]
      }
      
      // Generate the other entries.
      for transformID in blockStart..<blockEnd {
        // Load the column into the cache.
        var t = [Float](repeating: 0, count: blockSize)
        for rowID in blockStart..<transformID {
          let address = (rowID - blockStart) * blockSize + (transformID - blockStart)
          t[rowID - blockStart] = reflectorDotProducts[address]
        }
        
        // This GEMV operation could be transformed into a panel GEMM,
        // similarly to the method for accelerating classical Gram-Schmidt.
        var tt = [Float](repeating: 0, count: blockSize)
        let τ = currentMatrixT[transformID * problemSize + transformID]
        for rowID in blockStart..<blockEnd {
          // Multiply with the preceding submatrix.
          var dotProduct: Float = .zero
          for columnID in blockStart..<blockEnd {
            let matrixAddress = rowID * problemSize + columnID
            let matrixValue = currentMatrixT[matrixAddress]
            let vectorValue = t[columnID - blockStart]
            dotProduct += matrixValue * vectorValue
          }
          
          // Scale by the value on the diagonal.
          tt[rowID - blockStart] = -τ * dotProduct
        }
        
        // Store the column to main memory.
        for rowID in blockStart..<transformID {
          let address = rowID * problemSize + transformID
          currentMatrixT[address] = tt[rowID - blockStart]
        }
      }
      
      // Update by applying H**T to A(I:M,I+IB:N) from the right
      for columnID in 0..<problemSize {
        // V^H A
        var reflectorContributions = [Float](repeating: 0, count: blockSize)
        for transformID in blockStart..<blockEnd {
          var dotProduct: Float = .zero
          for rowID in 0..<problemSize {
            let matrixAddress = columnID * problemSize + rowID
            let matrixValue = matrix[matrixAddress]
            let vectorValue = panelReflectors[(transformID - blockStart) * problemSize + rowID]
            dotProduct += vectorValue * matrixValue
          }
          reflectorContributions[transformID - blockStart] = dotProduct
        }
        
        // T^H (V^H A)
        var TVA = [Float](repeating: 0, count: blockSize)
        for rowID in blockStart..<blockEnd {
          var dotProduct: Float = .zero
          for columnID in blockStart..<blockEnd {
            let addressT = columnID * problemSize + rowID
            let valueVA = reflectorContributions[columnID - blockStart]
            dotProduct += currentMatrixT[addressT] * valueVA
          }
          TVA[rowID - blockStart] = dotProduct
        }
        
        // V T^H V^H A
        for transformID in blockStart..<blockEnd {
          let dotProduct = TVA[transformID - blockStart]
          for rowID in 0..<problemSize {
            let matrixAddress = columnID * problemSize + rowID
            var matrixValue = matrix[matrixAddress]
            let vectorValue = panelReflectors[(transformID - blockStart) * problemSize + rowID]
            matrixValue -= dotProduct * vectorValue
            matrix[matrixAddress] = matrixValue
          }
        }
      }
      
      // Update by applying H**T to A(I:M,I+IB:N) from the left
      for columnID in 0..<problemSize {
        // V^H A
        var reflectorContributions = [Float](repeating: 0, count: blockSize)
        for transformID in blockStart..<blockEnd {
          var dotProduct: Float = .zero
          for rowID in 0..<problemSize {
            let matrixAddress = rowID * problemSize + columnID
            let matrixValue = matrix[matrixAddress]
            let vectorValue = panelReflectors[(transformID - blockStart) * problemSize + rowID]
            dotProduct += vectorValue * matrixValue
          }
          reflectorContributions[transformID - blockStart] = dotProduct
        }
        
        // T^H (V^H A)
        var TVA = [Float](repeating: 0, count: blockSize)
        for rowID in blockStart..<blockEnd {
          var dotProduct: Float = .zero
          for columnID in blockStart..<blockEnd {
            let addressT = columnID * problemSize + rowID
            let valueVA = reflectorContributions[columnID - blockStart]
            dotProduct += currentMatrixT[addressT] * valueVA
          }
          TVA[rowID - blockStart] = dotProduct
        }
        
        // V T^H V^H A
        for transformID in blockStart..<blockEnd {
          let dotProduct = TVA[transformID - blockStart]
          for rowID in 0..<problemSize {
            let matrixAddress = rowID * problemSize + columnID
            var matrixValue = matrix[matrixAddress]
            let vectorValue = panelReflectors[(transformID - blockStart) * problemSize + rowID]
            matrixValue -= dotProduct * vectorValue
            matrix[matrixAddress] = matrixValue
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
        currentTau[reflectorID] = panelTau[reflectorID - blockStart]
      }
    }
    
    return (currentReflectors, currentTau)
  }
}
