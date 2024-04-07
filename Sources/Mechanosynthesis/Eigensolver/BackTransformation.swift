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
    // let smallBlockSize = (blockSize + 1) / 2
    let smallBlockSize: Int = 2
    
    var rowOffset: Int = 1
    while rowOffset < problemSize {
      defer { rowOffset += blockSize }
      
      var blockStart = (problemSize - 1) / smallBlockSize * smallBlockSize
      while blockStart >= 0 {
        defer { blockStart -= smallBlockSize }
        
        // Establish bounds for 'rowID + elementID'.
        let remainingRowCount = max(0, problemSize - blockStart - rowOffset)
        let panelWidth = min(smallBlockSize, remainingRowCount)
        let panelHeight = min(blockSize + smallBlockSize, remainingRowCount)
        if panelHeight == 0 || panelWidth == 0 {
          continue
        }
        
        // Load the sweep into the cache.
        var reflectorBlock = [Float](
          repeating: .zero, 
          count: smallBlockSize * (blockSize + smallBlockSize))
        
        for sweepRelativeID in 0..<panelWidth {
          let sweepID = sweepRelativeID + blockStart
          let sweepMemoryOffset = sweepID * (problemSize + 1) + rowOffset
          let sweepCacheRow = smallBlockSize - 1 - sweepRelativeID
          var sweepCacheOffset = sweepCacheRow * (blockSize + smallBlockSize)
          sweepCacheOffset += sweepRelativeID
          
          let reflectorHeight = min(blockSize, remainingRowCount)
          for elementID in 0..<reflectorHeight {
            let matrixAddress = sweepMemoryOffset + elementID
            let matrixValue = bulgeReflectors[matrixAddress]
            reflectorBlock[sweepCacheOffset + elementID] = matrixValue
          }
        }
        
        // Create the T matrix using the 'WYTransform' API.
        var transformDesc = WYTransformDescriptor()
        transformDesc.dimension = SIMD2(blockSize + smallBlockSize, smallBlockSize)
        transformDesc.reflectorBlock = reflectorBlock
        let transform = WYTransform(descriptor: transformDesc)
        if blockSize == 4 {
//          print()
//          print("reflectors:", reflectorBlock)
//          print("tau", transform.tau)
//          print("blockStart:", blockStart)
//          print("rowOffset:", rowOffset)
//          print("problemSize:", problemSize)
//          
//          var reflectorDotProduct: Float = .zero
//          print("reflector dot product:", terminator: " ")
//          for elementID in 0..<blockSize + smallBlockSize {
//            let value0 = reflectorBlock[elementID]
//            let value1 = reflectorBlock[elementID + blockSize + smallBlockSize]
//            reflectorDotProduct += value0 * value1
//            print("(\(value0) * \(value1))", value0 * value1, terminator: " + ")
//          }
//          print("=", reflectorDotProduct)
        }
        
        #if true
        var VA = [Float](
          repeating: .zero, count: problemSize * smallBlockSize)
        
        for sweepRelativeID in 0..<smallBlockSize {
          let sweepBaseAddress = sweepRelativeID * (blockSize + smallBlockSize)
          for vectorID in 0..<problemSize {
            var vectorBaseAddress = vectorID * problemSize
            vectorBaseAddress += blockStart
            vectorBaseAddress += rowOffset
            
            var dotProduct: Float = .zero
            for elementID in 0..<panelHeight {
              let reflectorDatum = reflectorBlock[sweepBaseAddress + elementID]
              let vectorDatum = eigenvectors[vectorBaseAddress + elementID]
              dotProduct += reflectorDatum * vectorDatum
            }
            let dotProductAddress = vectorID * smallBlockSize + sweepRelativeID
            VA[dotProductAddress] = dotProduct
          }
        }
        
        let oldEigenvectors = eigenvectors
        if blockSize == 4 {
          //print("VA:", VA)
          //print("eigenvectors:", eigenvectors)
          
          
          
          //print("VA:", VA)
          //print("eigenvectors:", eigenvectors)
        }
        
        for vectorID in 0..<problemSize {
          let Tvalue = transform.tau[1]
          var VA0 = VA[vectorID * smallBlockSize]
          var VA1 = VA[vectorID * smallBlockSize + 1]
          //print(VA0, VA1, Tvalue, VA1 + VA0 * Tvalue)
          
          VA1 = VA1 + VA0 * Tvalue
          VA[vectorID * smallBlockSize + 1] = VA1
        }
        
        for sweepRelativeID in 0..<smallBlockSize {
          let sweepBaseAddress = sweepRelativeID * (blockSize + smallBlockSize)
          for vectorID in 0..<problemSize {
            var vectorBaseAddress = vectorID * problemSize
            vectorBaseAddress += blockStart
            vectorBaseAddress += rowOffset
            
            let dotProductAddress = vectorID * smallBlockSize + sweepRelativeID
            let dotProduct = VA[dotProductAddress]
            for elementID in 0..<panelHeight {
              let reflectorDatum = reflectorBlock[sweepBaseAddress + elementID]
              eigenvectors[vectorBaseAddress + elementID]
              -= reflectorDatum * dotProduct
            }
          }
        }
        
        let proposedNewEigenvectors = eigenvectors
        
        eigenvectors = oldEigenvectors
        VA = [Float](
          repeating: .zero, count: problemSize * smallBlockSize)
        
        // Back-transform the eigenvectors.
        for sweepRelativeID in 0..<smallBlockSize {
          let sweepBaseAddress = sweepRelativeID * (blockSize + smallBlockSize)
          for vectorID in 0..<problemSize {
            var vectorBaseAddress = vectorID * problemSize
            vectorBaseAddress += blockStart
            vectorBaseAddress += rowOffset
            
            var dotProduct: Float = .zero
            for elementID in 0..<panelHeight {
              let reflectorDatum = reflectorBlock[sweepBaseAddress + elementID]
              let vectorDatum = eigenvectors[vectorBaseAddress + elementID]
              dotProduct += reflectorDatum * vectorDatum
            }
            let dotProductAddress = vectorID * smallBlockSize + sweepRelativeID
            VA[dotProductAddress] = dotProduct
          }
          
          for vectorID in 0..<problemSize {
            var vectorBaseAddress = vectorID * problemSize
            vectorBaseAddress += blockStart
            vectorBaseAddress += rowOffset
            
            let dotProductAddress = vectorID * smallBlockSize + sweepRelativeID
            let dotProduct = VA[dotProductAddress]
            for elementID in 0..<panelHeight {
              let reflectorDatum = reflectorBlock[sweepBaseAddress + elementID]
              eigenvectors[vectorBaseAddress + elementID]
              -= reflectorDatum * dotProduct
            }
          }
        }
        
        if blockSize == 4 {
          for vectorID in 0..<problemSize {
            let Tvalue = transform.tau[1]
            let VA0 = VA[vectorID * smallBlockSize]
            let VA1 = VA[vectorID * smallBlockSize + 1]
            //print(VA0, VA1)
          }
          //print("VA:", VA)
          
//          print("old eigenvectors:         ", Array(oldEigenvectors[..<7]))
//          print("proposed new eigenvectors:", Array(proposedNewEigenvectors[..<7]))
//          print("eigenvectors:             ", Array(eigenvectors[..<7]))
          
          eigenvectors = proposedNewEigenvectors
        }
        
        
        #else
        // Apply the block reflector (A - V T^H V^H A).
        
        // TODO: Try operating on a single vector in isolation. It should
        // technically be possible to apply the block reflector via
        // matrix-vector multiplications.
        
        // TODO: Copy the vector's (used) data into a temporary cache. Apply 
        // the reflectors, then write back to memory. That may reveal an
        // incremental approach to integrating block reflections.
        
        // TODO: Perhaps debug by observing behavior of the second row. This
        // one has the simplest relationship between dot products with and w/o
        // the WY transform.
        // - set the small block size to 2
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
