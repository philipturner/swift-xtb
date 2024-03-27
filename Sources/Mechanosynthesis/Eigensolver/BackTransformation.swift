//
//  BackTransformation.swift
//
//
//  Created by Philip Turner on 3/11/24.
//

import Accelerate

extension Diagonalization {
  mutating func backTransform(
    bulgeSweeps: [BulgeSweep]
  ) {
    // If the matrix was already in tridiagonal form, there is no work to do.
    guard blockSize > 1 else {
      return
    }
    
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
        let sweep = bulgeSweeps[sweepID].data
        
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
      
      // Store the vector to main memory.
      for elementID in 0..<problemSize {
        let address = vectorID * problemSize + elementID
        eigenvectors[address] = vector[elementID]
      }
    }
  }
  
  mutating func backTransform(
    bandFormReflectors: [BandReflector]
  ) {
    for bandReflector in bandFormReflectors.reversed() {
      let panelReflectors = bandReflector.matrixV
      let T = bandReflector.matrixT
      
      // V^H A
      var VA = [Float](repeating: 0, count: problemSize * blockSize)
#if false
      for m in 0..<blockSize {
        for n in 0..<problemSize {
          var dotProduct: Float = .zero
          for k in 0..<problemSize {
            let lhsValue = panelReflectors[m * problemSize + k]
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
          &TRANSA, &TRANSB, &M, &N, &K, &ALPHA, panelReflectors, &LDA,
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
          &TRANSA, &TRANSB, &M, &N, &K, &ALPHA, panelReflectors, &LDA,
          TVA, &LDB, &BETA, &eigenvectors, &LDC)
      }
#endif
    }
  }
}
