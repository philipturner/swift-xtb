//
//  BandForm.swift
//
//
//  Created by Philip Turner on 3/11/24.
//

import Accelerate

extension Diagonalization {
  // Returns an array of reflector blocks.
  mutating func reduceToBandForm() -> [Float] {
    // Allocate a matrix to store the band reflectors.
    var bandReflectors = [Float](
      repeating: .zero, count: problemSize * problemSize)
    
    // Loop over the panels of the matrix.
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
      var reflectorBlock = [Float](
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
        
        // TODO: Refactor the application of reflectors to use BLAS and
        // recursive panel factorization.
        //
        // After this is finished, employ all possible constant-factor
        // improvements:
        // - exploiting symmetry,
        // - eliding multiplications by 0
        //
        // Precompute the WY transforms and check that there are no major
        // regressions. Parallelize everything as much as possible:
        // - parallelize the last two stages
        // - parallelize the bulge chasing
        //
        // The last few optimizations could be important in the future, but
        // perhaps not an economical use of time:
        // - pipeline the bulge chasing with panel factorization
        // - prepare the WY transforms simultaneously with other stages
        
        // Apply preceding reflectors (from this panel) to the column.
        for previousReflectorID in blockStart..<reflectorID {
          // Load the reflector into the cache.
          var reflector = [Float](repeating: 0, count: problemSize)
          for elementID in 0..<problemSize {
            let address = (
              previousReflectorID - blockStart) * problemSize + elementID
            reflector[elementID] = reflectorBlock[address]
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
        reflectorBlock.withContiguousMutableStorageIfAvailable { buffer in
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
      transformDesc.reflectorBlock = reflectorBlock
      let transform = WYTransform(descriptor: transformDesc)
      
      // Apply the reflector block from the left.
      do {
        // V^H A
        var VA = [Float](repeating: 0, count: problemSize * blockSize)
#if false
        for m in 0..<blockSize {
          for n in 0..<problemSize {
            var dotProduct: Float = .zero
            for k in 0..<problemSize {
              let lhsValue = reflectorBlock[m * problemSize + k]
              let rhsValue = matrix[k * problemSize + n]
              dotProduct += lhsValue * rhsValue
            }
            VA[n * blockSize + m] = dotProduct
          }
        }
#else
        var gemmDesc = GEMMDescriptor()
        gemmDesc.dimension = SIMD3(blockSize, problemSize, problemSize)
        reflectorBlock.withContiguousStorageIfAvailable {
          gemmDesc.leftOperand = $0.baseAddress!
          gemmDesc.leftOperandStride = problemSize
          gemmDesc.leftTransposeState = "T"
        }
        matrix.withContiguousStorageIfAvailable {
          gemmDesc.rightOperand = $0.baseAddress!
          gemmDesc.rightOperandStride = problemSize
        }
        VA.withContiguousMutableStorageIfAvailable {
          gemmDesc.accumulator = $0.baseAddress!
          gemmDesc.accumulatorStride = blockSize
        }
        GEMM(descriptor: gemmDesc)
        withExtendedLifetime(reflectorBlock) { }
        withExtendedLifetime(matrix) { }
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
        gemmDesc = GEMMDescriptor()
        gemmDesc.dimension = SIMD3(blockSize, problemSize, blockSize)
        transform.tau.withContiguousStorageIfAvailable {
          gemmDesc.leftOperand = $0.baseAddress!
          gemmDesc.leftOperandStride = blockSize
        }
        VA.withContiguousStorageIfAvailable {
          gemmDesc.rightOperand = $0.baseAddress!
          gemmDesc.rightOperandStride = blockSize
        }
        TVA.withContiguousMutableStorageIfAvailable {
          gemmDesc.accumulator = $0.baseAddress!
          gemmDesc.accumulatorStride = blockSize
        }
        GEMM(descriptor: gemmDesc)
        withExtendedLifetime(transform.tau) { }
        withExtendedLifetime(VA) { }
#endif
        
        // V (T^H V^H A)
#if false
        for m in 0..<problemSize {
          for n in 0..<problemSize {
            var dotProduct: Float = .zero
            for k in 0..<blockSize {
              let lhsValue = reflectorBlock[k * problemSize + m]
              let rhsValue = TVA[n * blockSize + k]
              dotProduct += lhsValue * rhsValue
            }
            matrix[m * problemSize + n] -= dotProduct
          }
        }
#else
        gemmDesc = GEMMDescriptor()
        gemmDesc.dimension = SIMD3(problemSize, problemSize, blockSize)
        gemmDesc.productScale = -1
        gemmDesc.accumulatorScale = 1
        TVA.withContiguousStorageIfAvailable {
          gemmDesc.leftOperand = $0.baseAddress!
          gemmDesc.leftOperandStride = blockSize
          gemmDesc.leftTransposeState = "T"
        }
        reflectorBlock.withContiguousStorageIfAvailable {
          gemmDesc.rightOperand = $0.baseAddress!
          gemmDesc.rightOperandStride = problemSize
          gemmDesc.rightTransposeState = "T"
        }
        matrix.withContiguousMutableStorageIfAvailable {
          gemmDesc.accumulator = $0.baseAddress!
          gemmDesc.accumulatorStride = problemSize
        }
        GEMM(descriptor: gemmDesc)
        withExtendedLifetime(TVA) { }
        withExtendedLifetime(reflectorBlock) { }
#endif
      }
      
      // Apply the reflector block from the right.
      do {
        // V^H A
        var VA = [Float](repeating: 0, count: problemSize * blockSize)
#if false
        for m in 0..<blockSize {
          for n in 0..<problemSize {
            var dotProduct: Float = .zero
            for k in 0..<problemSize {
              let lhsValue = reflectorBlock[m * problemSize + k]
              let rhsValue = matrix[n * problemSize + k]
              dotProduct += lhsValue * rhsValue
            }
            VA[n * blockSize + m] = dotProduct
          }
        }
#else
        var gemmDesc = GEMMDescriptor()
        gemmDesc.dimension = SIMD3(blockSize, problemSize, problemSize)
        reflectorBlock.withContiguousStorageIfAvailable {
          gemmDesc.leftOperand = $0.baseAddress!
          gemmDesc.leftOperandStride = problemSize
          gemmDesc.leftTransposeState = "T"
        }
        matrix.withContiguousStorageIfAvailable {
          gemmDesc.rightOperand = $0.baseAddress!
          gemmDesc.rightOperandStride = problemSize
        }
        VA.withContiguousMutableStorageIfAvailable {
          gemmDesc.accumulator = $0.baseAddress!
          gemmDesc.accumulatorStride = blockSize
        }
        GEMM(descriptor: gemmDesc)
        withExtendedLifetime(reflectorBlock) { }
        withExtendedLifetime(matrix) { }
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
        gemmDesc = GEMMDescriptor()
        gemmDesc.dimension = SIMD3(blockSize, problemSize, blockSize)
        transform.tau.withContiguousStorageIfAvailable {
          gemmDesc.leftOperand = $0.baseAddress!
          gemmDesc.leftOperandStride = blockSize
        }
        VA.withContiguousStorageIfAvailable {
          gemmDesc.rightOperand = $0.baseAddress!
          gemmDesc.rightOperandStride = blockSize
        }
        TVA.withContiguousMutableStorageIfAvailable {
          gemmDesc.accumulator = $0.baseAddress!
          gemmDesc.accumulatorStride = blockSize
        }
        GEMM(descriptor: gemmDesc)
        withExtendedLifetime(transform.tau) { }
        withExtendedLifetime(VA) { }
#endif
        
        // V (T^H V^H A)
#if false
        for m in 0..<problemSize {
          for n in 0..<problemSize {
            var dotProduct: Float = .zero
            for k in 0..<blockSize {
              let lhsValue = reflectorBlock[k * problemSize + m]
              let rhsValue = TVA[n * blockSize + k]
              dotProduct += lhsValue * rhsValue
            }
            matrix[n * problemSize + m] -= dotProduct
          }
        }
#else
        gemmDesc = GEMMDescriptor()
        gemmDesc.dimension = SIMD3(problemSize, problemSize, blockSize)
        gemmDesc.productScale = -1
        gemmDesc.accumulatorScale = 1
        reflectorBlock.withContiguousStorageIfAvailable {
          gemmDesc.leftOperand = $0.baseAddress!
          gemmDesc.leftOperandStride = problemSize
        }
        TVA.withContiguousStorageIfAvailable {
          gemmDesc.rightOperand = $0.baseAddress!
          gemmDesc.rightOperandStride = blockSize
        }
        matrix.withContiguousMutableStorageIfAvailable {
          gemmDesc.accumulator = $0.baseAddress!
          gemmDesc.accumulatorStride = problemSize
        }
        GEMM(descriptor: gemmDesc)
        withExtendedLifetime(TVA) { }
        withExtendedLifetime(reflectorBlock) { }
#endif
      }
      
      // Store the reflectors to main memory.
      for rowID in blockStart..<blockEnd {
        for columnID in 0..<problemSize {
          let matrixAddress = rowID * problemSize + columnID
          let panelAddress = (rowID - blockStart) * problemSize + columnID
          bandReflectors[matrixAddress] = reflectorBlock[panelAddress]
        }
      }
    }
    
    return bandReflectors
  }
}
