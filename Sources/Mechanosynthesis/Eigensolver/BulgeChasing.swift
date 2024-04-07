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
    
    // Loop over the bulge chasing sweeps.
    let sweepEnd = max(0, problemSize - 2)
    for sweepID in 0..<sweepEnd {
      var maxOperationID = (problemSize - 2) - (sweepID + 1)
      maxOperationID /= blockSize
      
      // Loop over the bulges within this sweep.
      for operationID in 0...maxOperationID {
        let startOfRow = operationID * blockSize
        let startOfColumn = sweepID * problemSize
        let startRowID = (sweepID + 1) + startOfRow
        
        // Create a reflector using the 'ReflectorGeneration' API.
        var generationDesc = ReflectorGenerationDescriptor()
        matrix.withContiguousStorageIfAvailable { buffer in
          let columnID = (sweepID + 1) + max(-1, startOfRow - blockSize)
          let matrixBaseAddress = columnID * problemSize + startRowID
          generationDesc.source = buffer.baseAddress! + matrixBaseAddress
        }
        
        // Find the address to begin writing data at.
        let reflectorBaseAddress = startOfColumn + (sweepID + 1) + startOfRow
        bulgeReflectors.withContiguousMutableStorageIfAvailable {
          generationDesc.destination = $0.baseAddress! + reflectorBaseAddress
        }
        
        // Determine the dimension of the reflector.
        let maxReflectorElementID = problemSize - sweepID - 1
        let endOfRow = min(startOfRow + blockSize, maxReflectorElementID)
        generationDesc.dimension = endOfRow - startOfRow
        ReflectorGeneration(descriptor: generationDesc)
        
        // Apply to the trailing submatrix.
        bulgeReflectors.withContiguousStorageIfAvailable {
          let reflector = $0.baseAddress! + reflectorBaseAddress
          let endRowID = min(startRowID + blockSize, problemSize)
          applyBulgeChase(
            reflector: reflector,
            startRowID: startRowID,
            endRowID: endRowID)
        }
      }
    }
    
    return bulgeReflectors
  }
  
  // Applies the reflector to the trailing submatrix.
  private mutating func applyBulgeChase(
    reflector: UnsafePointer<Float>,
    startRowID: Int,
    endRowID: Int
  ) {
    let rangeCount = endRowID - startRowID
    
    let startApplicationID: Int = max(startRowID - blockSize, 0)
    let endApplicationID: Int = min(endRowID + blockSize, problemSize)
    let dotProductCount = endApplicationID - startApplicationID
    var dotProducts = [Float](repeating: 0, count: dotProductCount)
    
    // Apply the reflector to the matrix, from the left.
    matrix.withContiguousMutableStorageIfAvailable {
      let matrixBaseAddress: Int =
      startApplicationID * problemSize + startRowID
      let matrix = $0.baseAddress! + matrixBaseAddress
      
      #if false
      for m in 0..<dotProductCount {
        for n in 0..<1 {
          var dotProduct: Float = .zero
          for k in 0..<rangeCount {
            let lhsValue = matrix[m * problemSize + k]
            let rhsValue = reflector[n * dotProductCount + k]
            dotProduct += lhsValue * rhsValue
          }
          dotProducts[m * 1 + n] = dotProduct
        }
      }
      #else
      do {
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
          &TRANSA, &TRANSB, &M, &N, &K, &ALPHA, matrix, &LDA,
          reflector, &LDB, &BETA, &dotProducts, &LDC)
      }
      #endif
      
      #if false
      for m in 0..<dotProductCount {
        for n in 0..<rangeCount {
          let lhsValue = dotProducts[m * 1]
          let rhsValue = reflector[n * 1]
          matrix[m * problemSize + n] -= lhsValue * rhsValue
        }
      }
      #else
      do {
        var M = Int32(truncatingIfNeeded: rangeCount)
        var N = Int32(truncatingIfNeeded: dotProductCount)
        var ALPHA = Float(-1)
        var INCX = Int32(1)
        var INCY = Int32(1)
        var LDA = Int32(truncatingIfNeeded: problemSize)
        sger_(
          &M, &N, &ALPHA, reflector, &INCX, dotProducts, &INCY,
          matrix, &LDA)
      }
      #endif
    }
    
    // Apply the reflector to the matrix, from the right.
    matrix.withContiguousMutableStorageIfAvailable {
      let matrixBaseAddress: Int = 
      startRowID * problemSize + startApplicationID
      let matrix = $0.baseAddress! + matrixBaseAddress
      
      #if false
      for m in 0..<dotProductCount {
        for n in 0..<1 {
          var dotProduct: Float = .zero
          for k in 0..<rangeCount {
            let lhsValue = matrix[k * problemSize + m]
            let rhsValue = reflector[n * dotProductCount + k]
            dotProduct += lhsValue * rhsValue
          }
          dotProducts[m * 1 + n] = dotProduct
        }
      }
      #else
      do {
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
          &TRANSA, &TRANSB, &M, &N, &K, &ALPHA, matrix, &LDA,
          reflector, &LDB, &BETA, &dotProducts, &LDC)
      }
      #endif
      
      #if false
      for m in 0..<rangeCount {
        for n in 0..<dotProductCount {
          let lhsValue = reflector[m]
          let rhsValue = dotProducts[n]
          matrix[m * problemSize + n] -= lhsValue * rhsValue
        }
      }
      #else
      do {
        var M = Int32(truncatingIfNeeded: dotProductCount)
        var N = Int32(truncatingIfNeeded: rangeCount)
        var ALPHA = Float(-1)
        var INCX = Int32(1)
        var INCY = Int32(1)
        var LDA = Int32(truncatingIfNeeded: problemSize)
        sger_(
          &M, &N, &ALPHA, dotProducts, &INCX, reflector, &INCY,
          matrix, &LDA)
      }
      #endif
    }
  }
}
