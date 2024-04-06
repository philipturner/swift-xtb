//
//  BackTransformation.swift
//
//
//  Created by Philip Turner on 3/11/24.
//

import Accelerate

extension Diagonalization {
  mutating func backTransform(
    bulgeChasingReflectors: [Float]
  ) {
    // Back-transform the eigenvectors.
    for vectorID in 0..<problemSize {
      // Load the vector into the cache.
      var vector = [Float](repeating: 0, count: problemSize)
      for elementID in 0..<problemSize {
        let address = vectorID * problemSize + elementID
        vector[elementID] = eigenvectors[address]
      }
      
      // TODO: Start optimizing this by breaking into sections modulo 32.
      // These will be pivot points, at which you reorder the reflector
      // elements.
      
      for sweepID in (0..<problemSize).reversed() {
        bulgeChasingReflectors.withContiguousStorageIfAvailable {
          let sweep = $0.baseAddress! + sweepID * problemSize
          var rowID: Int = sweepID + 1
          while rowID < problemSize {
            let nextRowID = min(rowID + blockSize, problemSize)
            let dotProductCount = nextRowID - rowID
            defer { rowID += blockSize }
            
            // Apply the reflector.
            var dotProduct: Float = .zero
            for elementID in 0..<dotProductCount {
              let reflectorDatum = sweep[rowID + elementID]
              let vectorDatum = vector[rowID + elementID]
              dotProduct += reflectorDatum * vectorDatum
            }
            for elementID in 0..<dotProductCount {
              let reflectorDatum = sweep[rowID + elementID]
              vector[rowID + elementID] -= reflectorDatum * dotProduct
            }
          }
        }
      }
      
      // Store the vector to main memory.
      for elementID in 0..<problemSize {
        let address = vectorID * problemSize + elementID
        eigenvectors[address] = vector[elementID]
      }
    }
  }
  
  mutating func backTransform(
    bandFormReflectors: [Float]
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
          reflectorBlock[panelAddress] = bandFormReflectors[matrixAddress]
        }
      }
      
      // Create the T matrix using the 'WYTransform' API.
      var transformDesc = WYTransformDescriptor()
      transformDesc = WYTransformDescriptor()
      transformDesc.dimension = SIMD2(problemSize, blockSize)
      transformDesc.reflectorBlock = reflectorBlock
      let transform = WYTransform(descriptor: transformDesc)
      
      // V^H A
      var VA = [Float](repeating: 0, count: problemSize * blockSize)
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
          &TRANSA, &TRANSB, &M, &N, &K, &ALPHA, reflectorBlock, &LDA,
          eigenvectors, &LDB, &BETA, &VA, &LDC)
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
            let lhsValue = reflectorBlock[k * problemSize + m]
            let rhsValue = TVA[n * blockSize + k]
            dotProduct += lhsValue * rhsValue
          }
          eigenvectors[n * problemSize + m] -= dotProduct
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
          &TRANSA, &TRANSB, &M, &N, &K, &ALPHA, reflectorBlock, &LDA,
          TVA, &LDB, &BETA, &eigenvectors, &LDC)
      }
#endif
    }
  }
}
