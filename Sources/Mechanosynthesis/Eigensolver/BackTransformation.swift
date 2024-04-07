//
//  BackTransformation.swift
//
//
//  Created by Philip Turner on 3/11/24.
//

import Accelerate

extension Diagonalization {
  mutating func backTransform(
    bulgeReflectors: [Float]
  ) {
#if false
    // Unoptimized implementation for debugging.
    for sweepID in (0..<problemSize).reversed() {
      bulgeReflectors.withContiguousStorageIfAvailable {
        let sweep = $0.baseAddress! + sweepID * problemSize
        
        var rowID: Int = sweepID + 1
        while rowID < problemSize {
          let nextRowID = min(rowID + blockSize, problemSize)
          let dotProductCount = nextRowID - rowID
          defer { rowID += blockSize }
          
          // Back-transform the eigenvectors.
          for vectorID in 0..<problemSize {
            let baseAddress = vectorID * problemSize
            
            var dotProduct: Float = .zero
            for elementID in 0..<dotProductCount {
              let reflectorDatum = sweep[rowID + elementID]
              let vectorDatum = eigenvectors[baseAddress + rowID + elementID]
              dotProduct += reflectorDatum * vectorDatum
            }
            for elementID in 0..<dotProductCount {
              let reflectorDatum = sweep[rowID + elementID]
              eigenvectors[baseAddress + rowID + elementID] -= reflectorDatum * dotProduct
            }
          }
        }
      }
    }
#else
    
    // The panels here are rectangular. The small block size is a heuristic to
    // minimize overhead, while keeping the growth in compute cost to <2x.
    let smallBlockSize = (blockSize + 3) / 4
    let smallProblemSize = blockSize + smallBlockSize
    
    // blockSize = 16, smallBlockSize = 8
    // 100 - 320
    // 250 - 2846
    // 500 - 16963
    // 750 - 49177
    // 1000 - 111923
    
    // blockSize = 16, smallBlockSize = 16
    // 100 - 201
    // 250 - 1247
    // 500 - 7789
    // 750 - 21661
    // 1000 - 49386
    
    // blockSize = 16, smallBlockSize = 32
    // 100 - 253
    // 250 - 1465
    // 500 - 7664
    // 750 - 21278
    // 1000 - 47899
    
    
    
    // blockSize = 32, smallBlockSize = 8
    // 100 - 174
    // 250 - 1613
    // 500 - 10143
    // 750 - 30517
    // 1000 - 69658
    
    // blockSize = 32, smallBlockSize = 16
    // 100 - 126
    // 250 - 779
    // 500 - 4907
    // 750 - 14310
    // 1000 - 33737
    
    // blockSize = 32, smallBlockSize = 32
    // 100 - 164
    // 250 - 754
    // 500 - 4091
    // 750 - 11673
    // 1000 - 28347
    
    // blockSize = 32, smallBlockSize = 64
    // 100 - 491
    // 250 - 1751
    // 500 - 7877
    // 750 - 19703
    // 1000 - 40990
    
    
    
    // blockSize = 64, smallBlockSize = 16
    // 100 - 90
    // 250 - 543
    // 500 - 3452
    // 750 - 10039
    // 1000 - 25686
    
    // blockSize = 64, smallBlockSize = 32
    // 100 - 109
    // 250 - 489
    // 500 - 3020
    // 750 - 7329
    // 1000 - 17091
    
    // blockSize = 64, smallBlockSize = 64
    // 100 - 236
    // 250 - 901
    // 500 - 4086
    // 750 - 10578
    // 1000 - 22205
    
    var rowOffset: Int = 1
    while rowOffset < problemSize {
      defer { rowOffset += blockSize }
      
      var blockStart = (problemSize - 1) / smallBlockSize * smallBlockSize
      while blockStart >= 0 {
        defer { blockStart -= smallBlockSize }
        
        // Establish bounds for 'rowID + elementID'.
        let remainingRowCount = max(0, problemSize - blockStart - rowOffset)
        let panelWidth = min(smallBlockSize, remainingRowCount)
        let panelHeight = min(smallProblemSize, remainingRowCount)
        if panelHeight == 0 || panelWidth == 0 {
          continue
        }
        
        // Load the sweep into the cache.
        var reflectorBlock = [Float](
          repeating: .zero,
          count: smallBlockSize * smallProblemSize)
        for m in 0..<panelWidth {
          let memorySweepID = m + blockStart
          let memoryBaseAddress = memorySweepID * (problemSize + 1) + rowOffset
          
          let cacheSweepID = smallBlockSize - 1 - m
          let cacheBaseAddress = cacheSweepID * smallProblemSize + m
          
          let reflectorLength = min(blockSize, remainingRowCount)
          for k in 0..<reflectorLength {
            let memoryValue = bulgeReflectors[memoryBaseAddress + k]
            reflectorBlock[cacheBaseAddress + k] = memoryValue
          }
        }
        
        // Create the T matrix using the 'WYTransform' API.
        var transformDesc = WYTransformDescriptor()
        transformDesc.dimension = SIMD2(smallProblemSize, smallBlockSize)
        transformDesc.reflectorBlock = reflectorBlock
        let transform = WYTransform(descriptor: transformDesc)
        
        // V^H A
        var VA = [Float](
          repeating: .zero, count: problemSize * smallBlockSize)
#if false
        for m in 0..<smallBlockSize {
          for n in 0..<problemSize {
            var dotProduct: Float = .zero
            for k in 0..<panelHeight {
              var vectorBaseAddress = n * problemSize
              vectorBaseAddress += blockStart
              vectorBaseAddress += rowOffset
              
              let lhsValue = reflectorBlock[m * smallProblemSize + k]
              let rhsValue = eigenvectors[vectorBaseAddress + k]
              dotProduct += lhsValue * rhsValue
            }
            VA[n * smallBlockSize + m] = dotProduct
          }
        }
#else
        var gemmDesc = GEMMDescriptor()
        gemmDesc.dimension = SIMD3(smallBlockSize, problemSize, panelHeight)
        reflectorBlock.withContiguousStorageIfAvailable {
          gemmDesc.leftOperand = $0.baseAddress!
          gemmDesc.leftOperandStride = smallProblemSize
          gemmDesc.leftTransposeState = "T"
        }
        eigenvectors.withContiguousStorageIfAvailable {
          gemmDesc.rightOperand = $0.baseAddress! + blockStart + rowOffset
          gemmDesc.rightOperandStride = problemSize
        }
        VA.withContiguousMutableStorageIfAvailable {
          gemmDesc.accumulator = $0.baseAddress!
          gemmDesc.accumulatorStride = smallBlockSize
        }
        GEMM(descriptor: gemmDesc)
        withExtendedLifetime(reflectorBlock) { }
        withExtendedLifetime(eigenvectors) { }
#endif
        
        // T^H (V^H A)
        var TVA = [Float](
          repeating: .zero, count: problemSize * smallBlockSize)
#if false
        for m in 0..<smallBlockSize {
          for n in 0..<problemSize {
            var dotProduct: Float = .zero
            for k in 0..<smallBlockSize {
              let lhsValue = transform.tau[k * smallBlockSize + m]
              let rhsValue = VA[n * smallBlockSize + k]
              dotProduct += lhsValue * rhsValue
            }
            TVA[n * smallBlockSize + m] = dotProduct
          }
        }
#else
        gemmDesc = GEMMDescriptor()
        gemmDesc.dimension = SIMD3(smallBlockSize, problemSize, smallBlockSize)
        transform.tau.withContiguousStorageIfAvailable {
          gemmDesc.leftOperand = $0.baseAddress!
          gemmDesc.leftOperandStride = smallBlockSize
        }
        VA.withContiguousStorageIfAvailable {
          gemmDesc.rightOperand = $0.baseAddress!
          gemmDesc.rightOperandStride = smallBlockSize
        }
        TVA.withContiguousMutableStorageIfAvailable {
          gemmDesc.accumulator = $0.baseAddress!
          gemmDesc.accumulatorStride = smallBlockSize
        }
        GEMM(descriptor: gemmDesc)
        withExtendedLifetime(transform.tau) { }
        withExtendedLifetime(VA) { }
#endif
        
        // V (T^H V^H A)
#if false
        for m in 0..<panelHeight {
          for n in 0..<problemSize {
            for k in 0..<smallBlockSize {
              var vectorBaseAddress = n * problemSize
              vectorBaseAddress += blockStart
              vectorBaseAddress += rowOffset
              
              let lhsValue = reflectorBlock[k * smallProblemSize + m]
              let rhsValue = TVA[n * smallBlockSize + k]
              eigenvectors[vectorBaseAddress + m] -= lhsValue * rhsValue
            }
          }
        }
#else
        gemmDesc = GEMMDescriptor()
        gemmDesc.dimension = SIMD3(panelHeight, problemSize, smallBlockSize)
        gemmDesc.productScale = -1
        gemmDesc.accumulatorScale = 1
        reflectorBlock.withContiguousStorageIfAvailable {
          gemmDesc.leftOperand = $0.baseAddress!
          gemmDesc.leftOperandStride = smallProblemSize
        }
        TVA.withContiguousStorageIfAvailable {
          gemmDesc.rightOperand = $0.baseAddress!
          gemmDesc.rightOperandStride = smallBlockSize
        }
        eigenvectors.withContiguousMutableStorageIfAvailable {
          gemmDesc.accumulator = $0.baseAddress! + blockStart + rowOffset
          gemmDesc.accumulatorStride = problemSize
        }
        GEMM(descriptor: gemmDesc)
        withExtendedLifetime(reflectorBlock) { }
        withExtendedLifetime(TVA) { }
#endif
      }
    }
#endif
  }
  
  mutating func backTransform(
    bandReflectors: [Float]
  ) {
    var blockStart = (problemSize - blockSize - 1) / blockSize * blockSize
    while blockStart >= 0 {
      // Adjust the loop end, to account for the factorization band offset.
      let blockEnd = min(blockStart + blockSize, problemSize - blockSize)
      defer { blockStart -= blockSize }
      
      // Reverse the order of the reflectors.
      var reflectorBlock = [Float](
        repeating: .zero, count: blockSize * problemSize)
      for rowID in blockStart..<blockEnd {
        var panelRowID = rowID - blockStart
        panelRowID = (blockSize - 1) - panelRowID
        
        for columnID in 0..<problemSize {
          let matrixAddress = rowID * problemSize + columnID
          let panelAddress = panelRowID * problemSize + columnID
          reflectorBlock[panelAddress] = bandReflectors[matrixAddress]
        }
      }
      
      // Create the T matrix using the 'WYTransform' API.
      var transformDesc = WYTransformDescriptor()
      transformDesc.dimension = SIMD2(problemSize, blockSize)
      transformDesc.reflectorBlock = reflectorBlock
      let transform = WYTransform(descriptor: transformDesc)
      
      // V^H A
      var VA = [Float](repeating: .zero, count: problemSize * blockSize)
#if false
      for m in 0..<blockSize {
        for n in 0..<problemSize {
          var dotProduct: Float = .zero
          for k in 0..<problemSize {
            let lhsValue = reflectorBlock[m * problemSize + k]
            let rhsValue = eigenvectors[n * problemSize + k]
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
      eigenvectors.withContiguousStorageIfAvailable {
        gemmDesc.rightOperand = $0.baseAddress!
        gemmDesc.rightOperandStride = problemSize
      }
      VA.withContiguousMutableStorageIfAvailable {
        gemmDesc.accumulator = $0.baseAddress!
        gemmDesc.accumulatorStride = blockSize
      }
      GEMM(descriptor: gemmDesc)
      withExtendedLifetime(reflectorBlock) { }
      withExtendedLifetime(eigenvectors) { }
#endif
      
      // T^H (V^H A)
      var TVA = [Float](repeating: .zero, count: problemSize * blockSize)
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
          eigenvectors[n * problemSize + m] -= dotProduct
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
      eigenvectors.withContiguousMutableStorageIfAvailable {
        gemmDesc.accumulator = $0.baseAddress!
        gemmDesc.accumulatorStride = problemSize
      }
      GEMM(descriptor: gemmDesc)
      withExtendedLifetime(reflectorBlock) { }
      withExtendedLifetime(TVA) { }
#endif
    }
  }
}
