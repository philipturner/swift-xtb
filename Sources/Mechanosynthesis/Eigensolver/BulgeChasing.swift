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
    var bulgeReflectorMatrix = [Float](
      repeating: 0, count: problemSize * problemSize)
    
    // Loop over the bulge chasing sweeps.
    let sweepEnd = max(0, problemSize - 2)
    for sweepID in 0..<sweepEnd {
      var maxOperationID = (problemSize - 2) - (sweepID + 1)
      maxOperationID /= blockSize
      
      // Loop over the bulges within this sweep.
      for operationID in 0...maxOperationID {
        let rowStart = operationID * blockSize
        var data: [Float]
        
        do {
          // The start and end vector ID must satisfy this condition:
          // endVectorID - startVectorID > 1
          let startColumnID = (sweepID + 1) + max(-1, rowStart - blockSize)
          let startVectorID = (sweepID + 1) + rowStart
          var endVectorID = (sweepID + 1) + (rowStart + blockSize)
          endVectorID = min(endVectorID, problemSize)
          
          // TODO: Extract some of the code from this function. Write the
          // bulge reflector directly into the matrix.
          data = applyBulgeChase(
            vectorID: startColumnID,
            startElementID: startVectorID,
            endElementID: endVectorID)
        }
        
        let maxReflectorElementID = problemSize - sweepID - 1
        let endReflectorElementID = min(
          rowStart + blockSize, maxReflectorElementID)
        let dotProductCount = endReflectorElementID - rowStart
        
        // Find the address to begin writing data at.
        let columnStart = sweepID * problemSize
        let reflectorBaseAddress = columnStart + (sweepID + 1) + rowStart
        for reflectorElementID in 0..<dotProductCount {
          let value = data[reflectorElementID]
          bulgeReflectorMatrix[reflectorBaseAddress + reflectorElementID] = value
        }
      }
    }
    
    return bulgeReflectorMatrix
  }
  
  // Returns a Householder reflector as an array allocation.
  private mutating func applyBulgeChase(
    vectorID: Int,
    startElementID: Int,
    endElementID: Int
  ) -> [Float] {
    let rangeCount = endElementID - startElementID
    
    // Load the row into the cache.
    var reflector = [Float](repeating: 0, count: rangeCount)
    for cacheElementID in 0..<rangeCount {
      let memoryElementID = startElementID + cacheElementID
      let matrixAddress = vectorID * problemSize + memoryElementID
      let matrixDatum = matrix[matrixAddress]
      reflector[cacheElementID] = matrixDatum
    }
    
    // Take the norm of the vector.
    var norm: Float = .zero
    for cacheElementID in 0..<rangeCount {
      let reflectorDatum = reflector[cacheElementID]
      norm += reflectorDatum * reflectorDatum
    }
    norm.formSquareRoot()
    
    // Predict the normalization factor.
    let oldSubdiagonal = reflector[0]
    let newSubdiagonal = norm * Float((oldSubdiagonal >= 0) ? -1 : 1)
    let tau = (newSubdiagonal - oldSubdiagonal) / newSubdiagonal
    
    // Check for NANs.
    var nanPresent = false
    let epsilon: Float = 2 * .leastNormalMagnitude
    if (newSubdiagonal - oldSubdiagonal).magnitude < epsilon {
      nanPresent = true
    }
    if newSubdiagonal.magnitude < epsilon {
      nanPresent = true
    }
    
    // Modify the vector, turning it into a reflector.
    for cacheElementID in 0..<rangeCount {
      var reflectorDatum = reflector[cacheElementID]
      if cacheElementID == 0 {
        reflectorDatum = 1
      } else {
        reflectorDatum /= oldSubdiagonal - newSubdiagonal
      }
      reflectorDatum *= tau.squareRoot()
      
      if nanPresent {
        reflectorDatum = 0
      }
      reflector[cacheElementID] = reflectorDatum
    }
    
    let startApplicationID: Int = max(startElementID - blockSize, 0)
    let endApplicationID: Int = min(endElementID + blockSize, problemSize)
    let dotProductCount = endApplicationID - startApplicationID
    var dotProducts = [Float](repeating: 0, count: dotProductCount)
    
    // Apply the reflector to the matrix, from the left.
    matrix.withContiguousMutableStorageIfAvailable {
      let matrixBaseAddress: Int =
      startApplicationID * problemSize + startElementID
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
        var TRANSA = CChar(Character("T").asciiValue!)
        var TRANSB = CChar(Character("N").asciiValue!)
        var M = Int32(dotProductCount)
        var N = Int32(1)
        var K = Int32(rangeCount)
        var ALPHA = Float(1)
        var LDA = Int32(problemSize)
        var BETA = Float(0)
        var LDB = Int32(dotProductCount)
        var LDC = Int32(dotProductCount)
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
        var M = Int32(rangeCount)
        var N = Int32(dotProductCount)
        var ALPHA = Float(-1)
        var INCX = Int32(1)
        var INCY = Int32(1)
        var LDA = Int32(problemSize)
        sger_(
          &M, &N, &ALPHA, reflector, &INCX, dotProducts, &INCY,
          matrix, &LDA)
      }
      #endif
    }
    
    // Apply the reflector to the matrix, from the right.
    matrix.withContiguousMutableStorageIfAvailable {
      let matrixBaseAddress: Int = 
      startElementID * problemSize + startApplicationID
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
        var TRANSA = CChar(Character("N").asciiValue!)
        var TRANSB = CChar(Character("N").asciiValue!)
        var M = Int32(dotProductCount)
        var N = Int32(1)
        var K = Int32(rangeCount)
        var ALPHA = Float(1)
        var LDA = Int32(problemSize)
        var BETA = Float(0)
        var LDB = Int32(dotProductCount)
        var LDC = Int32(dotProductCount)
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
        var M = Int32(dotProductCount)
        var N = Int32(rangeCount)
        var ALPHA = Float(-1)
        var INCX = Int32(1)
        var INCY = Int32(1)
        var LDA = Int32(problemSize)
        sger_(
          &M, &N, &ALPHA, dotProducts, &INCX, reflector, &INCY,
          matrix, &LDA)
      }
      #endif
    }
    
    // Store the reflector to main memory.
    return reflector
  }
}
