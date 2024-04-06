//
//  BandForm.swift
//
//
//  Created by Philip Turner on 3/11/24.
//

import Accelerate

struct BandReflector {
  // A matrix of reflectors, which describes a specific panel.
  var matrixV: [Float]
  
  // A matrix of coefficients to multiply dot products by.
  var matrixT: [Float]
}

extension Diagonalization {
  // Returns a matrix of reflectors.
  mutating func reduceToBandForm() -> [BandReflector] {
    var bandReflectors: [BandReflector] = []
    
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
        
        // TODO: Refactor the reflector application to use BLAS.
        
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
        
        // Create a reflector using the 'ReflectorGeneration' API.
        var generationDesc = ReflectorGenerationDescriptor()
        vector.withContiguousStorageIfAvailable { buffer in
          generationDesc.source = buffer.baseAddress! + bandOffset
        }
        panelReflectors.withContiguousMutableStorageIfAvailable { buffer in
          let offset = (reflectorID - blockStart) * problemSize + bandOffset
          generationDesc.destination = buffer.baseAddress! + offset
        }
        generationDesc.dimension = problemSize - bandOffset
        ReflectorGeneration(descriptor: generationDesc)
      }
      
      // Create the T matrix using the 'WYTransform' API.
      var transformDesc = WYTransformDescriptor()
      transformDesc = WYTransformDescriptor()
      transformDesc.dimension = SIMD2(problemSize, blockSize)
      transformDesc.reflectorBlock = panelReflectors
      var transform = WYTransform(descriptor: transformDesc)
      
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
              let lhsValue = transform.tau[k * blockSize + m]
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
            &TRANSA, &TRANSB, &M, &N, &K, &ALPHA, transform.tau, &LDA,
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
              let lhsValue = transform.tau[k * blockSize + m]
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
            &TRANSA, &TRANSB, &M, &N, &K, &ALPHA, transform.tau, &LDA,
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
      
      // TODO: Relocate this to the back-transformation stage. Make the T
      // matrices transient and store the V matrices in an n x n supermatrix.
      
      // Reverse the order of the reflectors.
      var reversedPanelReflectors = [Float](
        repeating: 0, count: blockSize * problemSize)
      for oldReflectorID in 0..<blockSize {
        let newReflectorID = blockSize - 1 - oldReflectorID
        for elementID in 0..<problemSize {
          let oldAddress = oldReflectorID * problemSize + elementID
          let newAddress = newReflectorID * problemSize + elementID
          reversedPanelReflectors[newAddress] = panelReflectors[oldAddress]
        }
      }
      
      // Create the T matrix using the 'WYTransform' API.
      transformDesc = WYTransformDescriptor()
      transformDesc.dimension = SIMD2(problemSize, blockSize)
      transformDesc.reflectorBlock = reversedPanelReflectors
      transform = WYTransform(descriptor: transformDesc)
      
      // Store the reflectors to main memory.
      bandReflectors.append(
        BandReflector(
          matrixV: reversedPanelReflectors,
          matrixT: transform.tau))
    }
    
    return bandReflectors
  }
}
