//
//  WYTransform.swift
//
//
//  Created by Philip Turner on 4/6/24.
//

import Accelerate

// A configuration for a compact WY transform.
struct WYTransformDescriptor {
  // The height of the panel (row count).
  var problemSize: Int?
  
  // The width of the panel (column count).
  var blockSize: Int?
  
  // The size of the internal blocking scheme.
  var smallBlockSize: Int?
  
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
    guard let problemSize = descriptor.problemSize,
          let blockSize = descriptor.blockSize,
          let smallBlockSize = descriptor.smallBlockSize,
          let reflectorBlock = descriptor.reflectorBlock else {
      fatalError("Descriptor not complete.")
    }
    
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
    
    // Initialize the T matrix to the identity matrix.
    tau = Array(repeating: .zero, count: blockSize * blockSize)
    for n in 0..<blockSize {
      tau[n * blockSize + n] = 1
    }
    
    // Fetch raw pointers for the matrices.
    let reflectorDotProductPointer = reflectorDotProducts
      .withContiguousStorageIfAvailable { $0.baseAddress! }!
    let tauPointer = tau
      .withContiguousMutableStorageIfAvailable { $0.baseAddress! }!
    
    // Fill in the T matrix in blocks.
    var blockStart: Int = .zero
    while blockStart < blockSize {
      let blockEnd = min(blockStart + smallBlockSize, blockSize)
      defer { blockStart += smallBlockSize }
      
      for k in blockStart..<blockEnd {
        let reflectorOffset = k * blockSize + (k + 1)
        let tauOffsetInput = k * blockSize
        let tauOffsetOutput = (k + 1) * blockSize
#if false
        for m in 0..<(blockEnd - k - 1) {
          for n in 0..<blockSize {
            let lhsValue = reflectorDotProducts[reflectorOffset + m]
            let rhsValue = tau[tauOffsetInput + n]
            tau[tauOffsetOutput + m * blockSize + n] -= lhsValue * rhsValue
          }
        }
#else
        let X = tauPointer + tauOffsetInput
        let Y = reflectorDotProductPointer + reflectorOffset
        let A = tauPointer + tauOffsetOutput
        
        var M = Int32(truncatingIfNeeded: blockSize)
        var N = Int32(truncatingIfNeeded: blockEnd - k - 1)
        var ALPHA = Float(-1)
        var INCX = Int32(1)
        var INCY = Int32(1)
        var LDA = Int32(truncatingIfNeeded: blockSize)
        sger_(
          &M,
          &N,
          &ALPHA,
          X, &INCX,
          UnsafeMutablePointer(mutating: Y), &INCY,
          A, &LDA)
#endif
      }
      
      guard blockSize > blockEnd else {
        continue
      }
#if false
      for m in 0..<blockSize {
        for n in 0..<blockSize - blockEnd {
          var dotProduct: Float = .zero
          for k in 0..<blockEnd - blockStart {
            var tauAddress = blockStart * blockSize
            tauAddress += k * blockSize + m
            
            var reflectorAddress = blockEnd * blockSize + blockStart
            reflectorAddress += n * blockSize + k
            
            let lhsValue = tau[tauAddress]
            let rhsValue = reflectorDotProducts[reflectorAddress]
            dotProduct += lhsValue * rhsValue
          }
          
          var tauAddress = blockEnd * blockSize
          tauAddress += n * blockSize + m
          tau[tauAddress] -= dotProduct
        }
      }
#else
      let dotProductOffset = blockEnd * blockSize + blockStart
      let A = UnsafePointer(tauPointer) + blockStart * blockSize
      let B = reflectorDotProductPointer + dotProductOffset
      let C = UnsafeMutablePointer(tauPointer) + blockEnd * blockSize
      
      var TRANSA = CChar(78) // N
      var TRANSB = CChar(78) // N
      var M = Int32(truncatingIfNeeded: blockSize)
      var N = Int32(truncatingIfNeeded: blockSize - blockEnd)
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
        UnsafeMutablePointer(mutating: A), &LDA,
        UnsafeMutablePointer(mutating: B), &LDB,
        &BETA,
        C, &LDC)
#endif
    }
    withExtendedLifetime(tau) { }
    withExtendedLifetime(reflectorDotProducts) { }
    
    transposeTau(blockSize: blockSize)
  }
  
  mutating func transposeTau(blockSize: Int) {
    var newTau = [Float](repeating: .zero, count: blockSize * blockSize)
    for rowID in 0..<blockSize {
      for columnID in 0..<blockSize {
        let oldValue = tau[rowID * blockSize + columnID]
        newTau[columnID * blockSize + rowID] = oldValue
      }
    }
    tau = newTau
  }
}
