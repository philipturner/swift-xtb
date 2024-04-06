//
//  GEMM.swift
//
//
//  Created by Philip Turner on 4/6/24.
//

import Accelerate

// A configuration for a general matrix-matrix multiplication.
struct GEMMDescriptor {
  // The M, N, and K arguments.
  var dimension: SIMD3<Int>?
  
  // The A, LDA, and TRANSA arguments.
  var leftOperand: UnsafePointer<Float>?
  var leftOperandStride: Int?
  var leftTransposeState: Character = "N"
  
  // The B, LDB, and TRANSB arguments.
  var rightOperand: UnsafePointer<Float>?
  var rightOperandStride: Int?
  var rightTransposeState: Character = "N"
  
  // The C and LDC arguments.
  var accumulator: UnsafeMutablePointer<Float>?
  var accumulatorStride: Int?
  
  // The 'alpha' and 'beta' arguments.
  var productScale: Float = 1
  var accumulatorScale: Float = 0
}

// Descriptor-based wrapper over `sgemm_` from Accelerate.
struct GEMM {
  // In typical API usage, one does not access the object's properties.
  @discardableResult
  init(descriptor: GEMMDescriptor) {
    guard let dimension = descriptor.dimension,
          let leftOperand = descriptor.leftOperand,
          let leftOperandStride = descriptor.leftOperandStride,
          let rightOperand = descriptor.rightOperand,
          let rightOperandStride = descriptor.rightOperandStride,
          let accumulator = descriptor.accumulator,
          let accumulatorStride = descriptor.accumulatorStride else {
      fatalError("Descriptor not complete.")
    }
    
    var TRANSA = CChar(descriptor.leftTransposeState.asciiValue!)
    var TRANSB = CChar(descriptor.rightTransposeState.asciiValue!)
    var M = Int32(dimension[0])
    var N = Int32(dimension[1])
    var K = Int32(dimension[2])
    var ALPHA = descriptor.productScale
    var LDA = Int32(leftOperandStride)
    var BETA = descriptor.accumulatorScale
    var LDB = Int32(rightOperandStride)
    var LDC = Int32(accumulatorStride)
    sgemm_(
      &TRANSA,
      &TRANSB,
      &M,
      &N,
      &K,
      &ALPHA,
      leftOperand, &LDA,
      rightOperand, &LDB,
      &BETA,
      accumulator, &LDC)
  }
}
