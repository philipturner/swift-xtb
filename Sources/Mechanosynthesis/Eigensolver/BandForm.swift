//
//  BandForm.swift
//
//
//  Created by Philip Turner on 3/11/24.
//

import Accelerate

extension Diagonalization {
  // TODO: Return an array of arrays, making it easier to pass the T matrix
  // from the reduction stage to the back-transformation stage.
  
  // Returns a matrix of reflectors.
  mutating func reduceToBandForm() -> (
    reflectors: [Float], tau: [Float]
  ) {
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
#if false
      for m in 0..<blockSize {
        for n in 0..<blockSize {
          var dotProduct: Float = .zero
          for k in 0..<problemSize {
            let lhsValue = panelReflectors[m * problemSize + k]
            let rhsValue = panelReflectors[n * problemSize + k]
            dotProduct += lhsValue * rhsValue
          }
          reflectorDotProducts[m * blockSize + n] = dotProduct
        }
      }
#else
      do {
        var TRANSA = CChar(Character("T").asciiValue!)
        var TRANSB = CChar(Character("N").asciiValue!)
        var M = Int32(blockSize)
        var N = Int32(blockSize)
        var K = Int32(problemSize)
        var ALPHA = Float(1)
        var LDA = Int32(problemSize)
        var BETA = Float(0)
        var LDB = Int32(problemSize)
        var LDC = Int32(blockSize)
        sgemm_(
          &TRANSA, &TRANSB, &M, &N, &K, &ALPHA, panelReflectors, &LDA,
          panelReflectors, &LDB, &BETA, &reflectorDotProducts, &LDC)
      }
#endif
      
      // Generate the T matrix. This could be saved and used to accelerate
      // the back-transformations.
      var T = [Float](repeating: 0, count: blockSize * blockSize)
      do {
        // Generate the diagonal entries.
        for n in 0..<blockSize {
          T[n * blockSize + n] = panelTau[n]
        }
        
        // Allocate cache memory for generating T.
        var tCache = [Float](repeating: 0, count: blockSize)
        var ttCache = [Float](repeating: 0, count: blockSize)
        
        // Generate the other entries.
        for n in 0..<blockSize {
          for m in 0..<blockSize {
            tCache[m] = reflectorDotProducts[m * blockSize + n]
          }
          
          // Multiply with the preceding submatrix.
          let τ = T[n * blockSize + n]
          for m in 0..<blockSize {
            var dotProduct: Float = .zero
            for k in 0..<blockSize {
              let matrixValue = T[m * blockSize + k]
              let vectorValue = tCache[k]
              dotProduct += matrixValue * vectorValue
            }
            ttCache[m] = -τ * dotProduct
          }
          
          for m in 0..<n {
            T[m * blockSize + n] = ttCache[m]
          }
        }
      }
      
      // MARK: - Update by applying H**T to A(I:M,I+IB:N) from the left
      
      do {
        // V^H A
        var VA = [Float](repeating: 0, count: problemSize * blockSize)
#if false
        for m in 0..<blockSize {
          for n in 0..<problemSize {
            var dotProduct: Float = .zero
            for k in 0..<problemSize {
              let lhsValue = panelReflectors[m * problemSize + k]
              let rhsValue = matrix[k * problemSize + n]
              dotProduct += lhsValue * rhsValue
            }
            VA[n * blockSize + m] = dotProduct
          }
        }
#else
        do {
          var TRANSA = CChar(Character("T").asciiValue!)
          var TRANSB = CChar(Character("T").asciiValue!)
          var M = Int32(blockSize)
          var N = Int32(problemSize)
          var K = Int32(problemSize)
          var ALPHA = Float(1)
          var LDA = Int32(problemSize)
          var BETA = Float(0)
          var LDB = Int32(problemSize)
          var LDC = Int32(blockSize)
          sgemm_(
            &TRANSA, &TRANSB, &M, &N, &K, &ALPHA, panelReflectors, &LDA,
            matrix, &LDB, &BETA, &VA, &LDC)
        }
#endif
        
        // T^H (V^H A)
        var TVA = [Float](repeating: 0, count: problemSize * blockSize)
#if false
        for m in 0..<blockSize {
          for n in 0..<problemSize {
            var dotProduct: Float = .zero
            for k in 0..<blockSize {
              let lhsValue = T[k * blockSize + m]
              let rhsValue = VA[n * blockSize + k]
              dotProduct += lhsValue * rhsValue
            }
            TVA[n * blockSize + m] = dotProduct
          }
        }
#else
        do {
          var TRANSA = CChar(Character("N").asciiValue!)
          var TRANSB = CChar(Character("N").asciiValue!)
          var M = Int32(blockSize)
          var N = Int32(problemSize)
          var K = Int32(blockSize)
          var ALPHA = Float(1)
          var LDA = Int32(blockSize)
          var BETA = Float(0)
          var LDB = Int32(blockSize)
          var LDC = Int32(blockSize)
          sgemm_(
            &TRANSA, &TRANSB, &M, &N, &K, &ALPHA, T, &LDA,
            VA, &LDB, &BETA, &TVA, &LDC)
        }
#endif
        
        // V (T^H V^H A)
        #if false
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
        #else
        do {
          var TRANSA = CChar(Character("T").asciiValue!)
          var TRANSB = CChar(Character("T").asciiValue!)
          var M = Int32(problemSize)
          var N = Int32(problemSize)
          var K = Int32(blockSize)
          var ALPHA = Float(-1)
          var LDA = Int32(blockSize)
          var BETA = Float(1)
          var LDB = Int32(problemSize)
          var LDC = Int32(problemSize)
          sgemm_(
            &TRANSA, &TRANSB, &M, &N, &K, &ALPHA, TVA, &LDA,
            panelReflectors, &LDB, &BETA, &matrix, &LDC)
        }
        #endif
      }
      
      // MARK: - Update by applying H**T to A(I:M,I+IB:N) from the right
      
      do {
        // V^H A
        var VA = [Float](repeating: 0, count: problemSize * blockSize)
#if false
        for m in 0..<blockSize {
          for n in 0..<problemSize {
            var dotProduct: Float = .zero
            for k in 0..<problemSize {
              let lhsValue = panelReflectors[m * problemSize + k]
              let rhsValue = matrix[n * problemSize + k]
              dotProduct += lhsValue * rhsValue
            }
            VA[n * blockSize + m] = dotProduct
          }
        }
#else
        do {
          var TRANSA = CChar(Character("T").asciiValue!)
          var TRANSB = CChar(Character("N").asciiValue!)
          var M = Int32(blockSize)
          var N = Int32(problemSize)
          var K = Int32(problemSize)
          var ALPHA = Float(1)
          var LDA = Int32(problemSize)
          var BETA = Float(0)
          var LDB = Int32(problemSize)
          var LDC = Int32(blockSize)
          sgemm_(
            &TRANSA, &TRANSB, &M, &N, &K, &ALPHA, panelReflectors, &LDA,
            matrix, &LDB, &BETA, &VA, &LDC)
        }
#endif
        
        // T^H (V^H A)
        var TVA = [Float](repeating: 0, count: problemSize * blockSize)
#if false
        for m in 0..<blockSize {
          for n in 0..<problemSize {
            var dotProduct: Float = .zero
            for k in 0..<blockSize {
              let lhsValue = T[k * blockSize + m]
              let rhsValue = VA[n * blockSize + k]
              dotProduct += lhsValue * rhsValue
            }
            TVA[n * blockSize + m] = dotProduct
          }
        }
#else
        do {
          var TRANSA = CChar(Character("N").asciiValue!)
          var TRANSB = CChar(Character("N").asciiValue!)
          var M = Int32(blockSize)
          var N = Int32(problemSize)
          var K = Int32(blockSize)
          var ALPHA = Float(1)
          var LDA = Int32(blockSize)
          var BETA = Float(0)
          var LDB = Int32(blockSize)
          var LDC = Int32(blockSize)
          sgemm_(
            &TRANSA, &TRANSB, &M, &N, &K, &ALPHA, T, &LDA,
            VA, &LDB, &BETA, &TVA, &LDC)
        }
#endif
        
        // V (T^H V^H A)
#if false
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
#else
        do {
          var TRANSA = CChar(Character("N").asciiValue!)
          var TRANSB = CChar(Character("N").asciiValue!)
          var M = Int32(problemSize)
          var N = Int32(problemSize)
          var K = Int32(blockSize)
          var ALPHA = Float(-1)
          var LDA = Int32(problemSize)
          var BETA = Float(1)
          var LDB = Int32(blockSize)
          var LDC = Int32(problemSize)
          sgemm_(
            &TRANSA, &TRANSB, &M, &N, &K, &ALPHA, panelReflectors, &LDA,
            TVA, &LDB, &BETA, &matrix, &LDC)
        }
#endif
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
