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
      
      // TODO: Reformulate every statement below as a matrix multiplication in
      // scalar code. After doing that, replace each matrix multiplication with
      // a call to Accelerate.
      
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
      for transformID in 0..<blockEnd - blockStart {
        let diagonalAddress = (blockStart + transformID) * problemSize + (blockStart + transformID)
        let dotProductAddress = transformID * blockSize + transformID
        currentMatrixT[diagonalAddress] = panelTau[transformID]
      }
      
      // Generate the other entries.
      for transformID in 0..<blockEnd - blockStart {
        // Load the column into the cache.
        var t = [Float](repeating: 0, count: blockSize)
        for rowID in 0..<transformID {
          let address = rowID * blockSize + transformID
          t[rowID] = reflectorDotProducts[address]
        }
        
        // This GEMV operation could be transformed into a panel GEMM,
        // similarly to the method for accelerating classical Gram-Schmidt.
        var tt = [Float](repeating: 0, count: blockSize)
        let τ = currentMatrixT[(blockStart + transformID) * problemSize + (blockStart + transformID)]
        for rowID in 0..<blockEnd - blockStart {
          // Multiply with the preceding submatrix.
          var dotProduct: Float = .zero
          for columnID in 0..<blockEnd - blockStart {
            let matrixAddress = (blockStart + rowID) * problemSize + (blockStart + columnID)
            let matrixValue = currentMatrixT[matrixAddress]
            let vectorValue = t[columnID]
            dotProduct += matrixValue * vectorValue
          }
          
          // Scale by the value on the diagonal.
          tt[rowID] = -τ * dotProduct
        }
        
        // Store the column to main memory.
        for rowID in 0..<transformID {
          let address = (blockStart + rowID) * problemSize + (blockStart + transformID)
          currentMatrixT[address] = tt[rowID]
        }
      }
      
      // Use this to incrementally debug the removal of the excessive memory
      // allocation for T.
      var compressedT = [Float](repeating: 0, count: blockSize * blockSize)
      for rowID in 0..<blockEnd - blockStart {
        for columnID in 0..<blockEnd - blockStart {
          compressedT[rowID * blockSize + columnID] = currentMatrixT[(blockStart + rowID) * problemSize + (blockStart + columnID)]
        }
      }
      
      // MARK: - Update by applying H**T to A(I:M,I+IB:N) from the right
      
      do {
        // V^H A
        var VA = [Float](repeating: 0, count: problemSize * blockSize)
        for m in 0..<problemSize {
          for n in 0..<blockSize {
            var dotProduct: Float = .zero
            for k in 0..<problemSize {
              let lhsValue = panelReflectors[n * problemSize + k]
              let rhsValue = matrix[m * problemSize + k]
              dotProduct += lhsValue * rhsValue
            }
            VA[m * blockSize + n] = dotProduct
          }
        }
        
        // T^H (V^H A)
        var TVA = [Float](repeating: 0, count: problemSize * blockSize)
        for m in 0..<problemSize {
          for n in 0..<blockSize {
            var dotProduct: Float = .zero
            for k in 0..<blockSize {
              let lhsValue = compressedT[k * blockSize + n]
              let rhsValue = VA[m * blockSize + k]
              dotProduct += lhsValue * rhsValue
            }
            TVA[m * blockSize + n] = dotProduct
          }
        }
        
        // V (T^H V^H A)
        for m in 0..<problemSize {
          for n in 0..<problemSize {
            var dotProduct: Float = .zero
            for k in 0..<blockSize {
              let lhsValue = panelReflectors[k * problemSize + m]
              let rhsValue = TVA[n * blockSize + k]
              dotProduct += lhsValue * rhsValue
            }
            matrix[n * problemSize + m] -= dotProduct
          }
        }
      }
      
      // MARK: - Update by applying H**T to A(I:M,I+IB:N) from the left
      
      do {
        // V^H A
        var VA = [Float](repeating: 0, count: problemSize * blockSize)
        for m in 0..<problemSize {
          for n in 0..<blockSize {
            var dotProduct: Float = .zero
            for k in 0..<problemSize {
              let lhsValue = panelReflectors[n * problemSize + k]
              let rhsValue = matrix[k * problemSize + m]
              dotProduct += lhsValue * rhsValue
            }
            VA[m * blockSize + n] = dotProduct
          }
        }
        
        // T^H (V^H A)
        var TVA = [Float](repeating: 0, count: problemSize * blockSize)
        for m in 0..<problemSize {
          for n in 0..<blockSize {
            var dotProduct: Float = .zero
            for k in 0..<blockSize {
              let lhsValue = compressedT[k * blockSize + n]
              let rhsValue = VA[m * blockSize + k]
              dotProduct += lhsValue * rhsValue
            }
            TVA[m * blockSize + n] = dotProduct
          }
        }
        
        // V (T^H V^H A)
        for m in 0..<problemSize {
          for n in 0..<problemSize {
            var dotProduct: Float = .zero
            for k in 0..<blockSize {
              let lhsValue = panelReflectors[k * problemSize + m]
              let rhsValue = TVA[n * blockSize + k]
              dotProduct += lhsValue * rhsValue
            }
            matrix[m * problemSize + n] -= dotProduct
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
