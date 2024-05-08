//
//  BulgeChasing.swift
//
//
//  Created by Philip Turner on 3/11/24.
//

import Accelerate

extension Diagonalization {
  // Returns a sequence of reflectors.
  mutating func chaseBulges() -> [Float] {
    // Allocate a matrix to store the bulge reflectors.
    var bulgeReflectors = [Float](
      repeating: .zero, count: problemSize * problemSize)
    let reflectors = bulgeReflectors
      .withContiguousMutableStorageIfAvailable { $0.baseAddress! }!
    let dotProducts: UnsafeMutablePointer<Float> =
      .allocate(capacity: 3 * blockSize)
    defer { dotProducts.deallocate() }
    
    // Loop over the bulge chasing sweeps.
    let sweepEnd = max(0, problemSize - 2)
    for sweepID in 0..<sweepEnd {
      var maxOperationID = (problemSize - 2) - (sweepID + 1)
      maxOperationID /= blockSize
      
      // Loop over the bulges within this sweep.
      for operationID in 0...maxOperationID {
        let startOfRow = operationID &* blockSize
        let startOfColumn = sweepID &* problemSize
        let startRowID = (sweepID + 1) &+ startOfRow
        
        // Create a reflector using the 'ReflectorGeneration' API.
        var generationDesc = ReflectorGenerationDescriptor()
        let columnID = (sweepID + 1) &+ max(-1, startOfRow &- blockSize)
        let matrixBaseAddress = columnID &* problemSize &+ startRowID
        let matrixSource = UnsafePointer(matrixPointer! + matrixBaseAddress)
        generationDesc.source = matrixSource
        
        // Find the address to begin writing data at.
        let reflectorBaseAddress = startOfColumn &+ (sweepID + 1) &+ startOfRow
        let reflector = reflectors + reflectorBaseAddress
        generationDesc.destination = reflector
        
        // Determine the dimension of the reflector.
        let maxReflectorElementID = problemSize &- (sweepID + 1)
        let endOfRow = min(startOfRow &+ blockSize, maxReflectorElementID)
        generationDesc.dimension = endOfRow &- startOfRow
        ReflectorGeneration(descriptor: generationDesc)
        
        // Apply to the trailing submatrix.
        let endRowID = min(startRowID &+ blockSize, problemSize)
        applyBulgeChase(
          reflector: reflector,
          dotProducts: dotProducts,
          startRowID: startRowID,
          endRowID: endRowID)
      }
    }
    
    return bulgeReflectors
  }
  
  // Applies the reflector to the trailing submatrix.
  @_transparent
  private mutating func applyBulgeChase(
    reflector: UnsafePointer<Float>,
    dotProducts: UnsafeMutablePointer<Float>,
    startRowID: Int,
    endRowID: Int
  ) {
    let startApplicationID = max(startRowID &- blockSize, 0)
    let endApplicationID = min(endRowID &+ blockSize, problemSize)
    let dotProductCount = endApplicationID &- startApplicationID
    let rangeCount = endRowID &- startRowID
    
    // Apply the reflector to the matrix, from the left.
    do {
      let matrixOffset = startApplicationID &* problemSize &+ startRowID
      let A = matrixPointer! + matrixOffset
      let B = reflector
      let C = dotProducts
#if false
      for m in 0..<dotProductCount {
        for n in 0..<1 {
          var dotProduct: Float = .zero
          for k in 0..<rangeCount {
            let lhsValue = matrix[matrixOffset + m * problemSize + k]
            let rhsValue = reflector[n * dotProductCount + k]
            dotProduct += lhsValue * rhsValue
          }
          dotProducts[m * 1 + n] = dotProduct
        }
      }
#else
      var TRANSA = CChar(84) // T
      var TRANSB = CChar(78) // N
      var M = Int32(truncatingIfNeeded: dotProductCount)
      var N = Int32(1)
      var K = Int32(truncatingIfNeeded: rangeCount)
      var ALPHA = Float(1)
      var LDA = Int32(truncatingIfNeeded: problemSize)
      var BETA = Float(0)
      var LDB = Int32(truncatingIfNeeded: dotProductCount)
      var LDC = Int32(truncatingIfNeeded: dotProductCount)
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
    
    do {
      let matrixOffset = startApplicationID &* problemSize &+ startRowID
      let X = reflector
      let Y = dotProducts
      let A = matrixPointer! + matrixOffset
#if false
      for m in 0..<dotProductCount {
        for n in 0..<rangeCount {
          let lhsValue = dotProducts[m * 1]
          let rhsValue = reflector[n * 1]
          matrix[matrixOffset + m * problemSize + n] -= lhsValue * rhsValue
        }
      }
#else
      var M = Int32(truncatingIfNeeded: rangeCount)
      var N = Int32(truncatingIfNeeded: dotProductCount)
      var ALPHA = Float(-1)
      var INCX = Int32(1)
      var INCY = Int32(1)
      var LDA = Int32(truncatingIfNeeded: problemSize)
      sger_(
        &M,
        &N,
        &ALPHA,
        X, &INCX,
        Y, &INCY,
        A, &LDA)
#endif
    }
    
    // Apply the reflector to the matrix, from the right.
    do {
      let matrixOffset = startRowID &* problemSize &+ startApplicationID
      let A = matrixPointer! + matrixOffset
      let B = reflector
      let C = dotProducts
#if false
      for m in 0..<dotProductCount {
        for n in 0..<1 {
          var dotProduct: Float = .zero
          for k in 0..<rangeCount {
            let lhsValue = matrix[matrixOffset + k * problemSize + m]
            let rhsValue = reflector[n * dotProductCount + k]
            dotProduct += lhsValue * rhsValue
          }
          dotProducts[m * 1 + n] = dotProduct
        }
      }
#else
      var TRANSA = CChar(78) // N
      var TRANSB = CChar(78) // N
      var M = Int32(truncatingIfNeeded: dotProductCount)
      var N = Int32(1)
      var K = Int32(truncatingIfNeeded: rangeCount)
      var ALPHA = Float(1)
      var LDA = Int32(truncatingIfNeeded: problemSize)
      var BETA = Float(0)
      var LDB = Int32(truncatingIfNeeded: dotProductCount)
      var LDC = Int32(truncatingIfNeeded: dotProductCount)
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
    
    do {
      let matrixOffset = startRowID &* problemSize &+ startApplicationID
      let X = dotProducts
      let Y = reflector
      let A = matrixPointer! + matrixOffset
#if false
      for m in 0..<rangeCount {
        for n in 0..<dotProductCount {
          let lhsValue = reflector[m]
          let rhsValue = dotProducts[n]
          matrix[matrixOffset + m * problemSize + n] -= lhsValue * rhsValue
        }
      }
#else
      var M = Int32(truncatingIfNeeded: dotProductCount)
      var N = Int32(truncatingIfNeeded: rangeCount)
      var ALPHA = Float(-1)
      var INCX = Int32(1)
      var INCY = Int32(1)
      var LDA = Int32(truncatingIfNeeded: problemSize)
      sger_(
        &M,
        &N,
        &ALPHA,
        X, &INCX,
        Y, &INCY,
        A, &LDA)
#endif
    }
  }
}
