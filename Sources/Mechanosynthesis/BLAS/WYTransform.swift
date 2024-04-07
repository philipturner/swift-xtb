//
//  WYTransform.swift
//
//
//  Created by Philip Turner on 4/6/24.
//

import Accelerate

// A configuration for a compact WY transform.
struct WYTransformDescriptor {
  // The dimensions (rows, columns) of the block. The first entry should be the
  // leading dimension.
  var dimension: SIMD2<Int>?
  
  // The block of reflectors, stored in column-major order.
  var reflectorBlock: [Float]?
  
  @_transparent
  init() { }
}

// Applies the 'compact WY transform' to a block of reflectors.
struct WYTransform {
  // The T matrix.
  // - dimensions: dimension[1] x dimension[1]
  // - order: unknown
  var tau: [Float]
  
  @_transparent
  init(descriptor: WYTransformDescriptor) {
    guard let dimension = descriptor.dimension,
          let reflectorBlock = descriptor.reflectorBlock else {
      fatalError("Descriptor not complete.")
    }
    let problemSize = dimension[0]
    let blockSize = dimension[1]
    
    // Generate an overlap matrix for the reflectors.
    var reflectorDotProducts = [Float](
      repeating: .zero, count: blockSize * blockSize)
#if false
    for m in 0..<blockSize {
      for n in 0..<blockSize {
        var dotProduct: Float = .zero
        for k in 0..<problemSize {
          let lhsValue = reflectorBlock[m * problemSize + k]
          let rhsValue = reflectorBlock[n * problemSize + k]
          dotProduct += lhsValue * rhsValue
        }
        reflectorDotProducts[m * blockSize + n] = dotProduct
      }
    }
#else
    var gemmDesc = GEMMDescriptor()
    gemmDesc.dimension = SIMD3(blockSize, blockSize, problemSize)
    reflectorBlock.withContiguousStorageIfAvailable {
      gemmDesc.leftOperand = $0.baseAddress!
      gemmDesc.leftOperandStride = problemSize
      gemmDesc.leftTransposeState = "T"
    }
    reflectorBlock.withContiguousStorageIfAvailable {
      gemmDesc.rightOperand = $0.baseAddress!
      gemmDesc.rightOperandStride = problemSize
    }
    reflectorDotProducts.withContiguousMutableStorageIfAvailable {
      gemmDesc.accumulator = $0.baseAddress!
      gemmDesc.accumulatorStride = blockSize
    }
    GEMM(descriptor: gemmDesc)
#endif
    withExtendedLifetime(reflectorBlock) { }
    
    // Use a heuristic for the recursive block size.
    let smallBlockSize = max(1, min((blockSize + 3) / 4, blockSize))
    
    // Make the T matrix technically initialized, according to the Swift
    // compiler.
    tau = []
    for _ in 0..<4 {
      // Initialize the T matrix to the identity matrix.
      tau = Array(repeating: .zero, count: blockSize * blockSize)
      for n in 0..<blockSize {
        tau[n * blockSize + n] = 1
      }
      
      let reflectorDotProductPointer = reflectorDotProducts
        .withContiguousStorageIfAvailable { $0.baseAddress! }!
      let tauPointer = tau
        .withContiguousMutableStorageIfAvailable { $0.baseAddress! }!
      
      // Fill in the T matrix in blocks.
      var blockStart: Int = .zero
      while blockStart < blockSize {
        let blockEnd = min(blockStart + smallBlockSize, blockSize)
        defer { blockStart += smallBlockSize }
        
        for m in 0..<blockSize {
          for k in blockStart..<blockEnd {
            for n in (k + 1)..<blockEnd {
              let matrixValue = tau[m * blockSize + k]
              let vectorValue = reflectorDotProducts[k * blockSize + n]
              tau[m * blockSize + n] -= matrixValue * vectorValue
            }
          }
        }
        
#if false
        for m in 0..<blockSize - blockEnd {
          for n in 0..<blockSize {
            var dotProduct: Float = .zero
            for k in 0..<blockEnd - blockStart {
              var reflectorAddress = blockEnd * blockSize + blockStart
              reflectorAddress += m * blockSize + k
              
              var tauAddress = blockStart
              tauAddress += n * blockSize + k
              
              let lhsValue = reflectorDotProducts[reflectorAddress]
              let rhsValue = tau[tauAddress]
              dotProduct += lhsValue * rhsValue
            }
            
            var tauAddress = blockEnd
            tauAddress += n * blockSize + m
            tau[tauAddress] -= dotProduct
          }
        }
#else
        let dotProductOffset = blockEnd * blockSize + blockStart
        let A = reflectorDotProductPointer + dotProductOffset
        let B = UnsafePointer(tauPointer) + blockStart
        let C = UnsafeMutablePointer(tauPointer) + blockEnd
        
        var TRANSA = CChar(84) // T
        var TRANSB = CChar(78) // N
        var M = Int32(truncatingIfNeeded: blockSize - blockEnd)
        var N = Int32(truncatingIfNeeded: blockSize)
        var K = Int32(truncatingIfNeeded: blockEnd - blockStart)
        var ALPHA = Float(-1)
        var LDA = Int32(truncatingIfNeeded: blockSize)
        var BETA = Float(1)
        var LDB = Int32(truncatingIfNeeded: blockSize)
        var LDC = Int32(truncatingIfNeeded: blockSize)
        sgemm_(
          &TRANSA,
          &TRANSB,
          &M,
          &N,
          &K,
          &ALPHA,
          A, &LDA,
          B, &LDB,
          &BETA,
          C, &LDC)
#endif
      }
      
      withExtendedLifetime(tau) { }
      withExtendedLifetime(reflectorDotProducts) { }
    }
  }
}
