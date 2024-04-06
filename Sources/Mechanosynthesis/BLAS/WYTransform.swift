//
//  WYTransform.swift
//
//
//  Created by Philip Turner on 4/6/24.
//

// A configuration for a compact WY transform.
struct WYTransformDescriptor {
  // The dimensions (rows, columns) of the block. The first entry should be the
  // leading dimension.
  var dimension: SIMD2<Int>?
  
  // The block of reflectors, stored in column-major order.
  var reflectorBlock: [Float]?
  
  // The scalar 'tau' values for each reflector.
  var tauBlock: [Float]?
}

// Applies the _compact WY transform_ to a block of reflectors.
struct WYTransform {
  // The T matrix.
  // - dimensions: dimension[1] x dimension[1]
  // - layout: unknown
  var tau: [Float]
  
  init(descriptor: WYTransformDescriptor) {
    guard let dimension = descriptor.dimension,
          let reflectorBlock = descriptor.reflectorBlock,
          let tauBlock = descriptor.tauBlock else {
      fatalError("Descriptor not complete.")
    }
    let problemSize = dimension[0]
    let blockSize = dimension[1]
    tau = Array(repeating: .zero, count: blockSize * blockSize)
    
    // Generate an overlap matrix for the reflectors.
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
    
    // Generate the diagonal entries.
    for n in 0..<blockSize {
      tau[n * blockSize + n] = tauBlock[n]
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
      let τ = tau[n * blockSize + n]
      for m in 0..<blockSize {
        var dotProduct: Float = .zero
        for k in 0..<blockSize {
          let matrixValue = tau[m * blockSize + k]
          let vectorValue = tCache[k]
          dotProduct += matrixValue * vectorValue
        }
        ttCache[m] = -τ * dotProduct
      }
      
      for m in 0..<n {
        tau[m * blockSize + n] = ttCache[m]
      }
    }
  }
}
