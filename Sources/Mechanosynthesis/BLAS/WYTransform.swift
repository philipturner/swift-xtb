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
    tau = Array(repeating: .zero, count: blockSize * blockSize)
    
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
    for n in 0..<blockSize {
      tau[n * blockSize + n] = 1
    }
    
    // Generate the other entries.
    for n in 0..<blockSize {
      for m in 0..<n {
        var dotProduct: Float = .zero
        for k in 0..<blockSize {
          let matrixValue = tau[m * blockSize + k]
          let vectorValue = reflectorDotProducts[k * blockSize + n]
          dotProduct += matrixValue * vectorValue
        }
        tau[m * blockSize + n] = -dotProduct
      }
    }
  }
}
