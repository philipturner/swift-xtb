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
        
        var reflector: [Float]
        do {
          // TODO: Change 'vectorID' and 'elementID' to 'rowID'.
          
          // The start and end vector ID must satisfy this condition:
          // endVectorID - startVectorID > 1
          let startColumnID = (sweepID + 1) + max(-1, rowStart - blockSize)
          let startVectorID = (sweepID + 1) + rowStart
          let endVectorID = min(startVectorID + blockSize, problemSize)
          
          // TODO: Refactor the code in this function. Write the bulge
          // reflector directly into the matrix.
          reflector = createBulgeReflector(
            matrixBaseAddress: startColumnID * problemSize + startVectorID,
            startElementID: startVectorID,
            endElementID: endVectorID)
          
          // TODO: Refactor the code in this function. Read the bulge
          // reflector directly from the matrix.
          applyBulgeChase(
            reflector: reflector,
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
          let value = reflector[reflectorElementID]
          bulgeReflectorMatrix[reflectorBaseAddress + reflectorElementID] = value
        }
      }
    }
    
    return bulgeReflectorMatrix
  }
  
  // Returns a Householder reflector as an array allocation.
  private func createBulgeReflector(
    matrixBaseAddress: Int,
    startElementID: Int,
    endElementID: Int
  ) -> [Float] {
    let rangeCount = endElementID - startElementID
    
    // Load the row into the cache.
    var reflector = [Float](repeating: 0, count: blockSize)
    for elementID in 0..<rangeCount {
      let matrixAddress = matrixBaseAddress + elementID
      let matrixDatum = matrix[matrixAddress]
      reflector[elementID] = matrixDatum
    }
    
    // Take the norm of the vector.
    var norm: Float = .zero
    for elementID in 0..<blockSize {
      let reflectorDatum = reflector[elementID]
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
    for elementID in 0..<blockSize {
      var reflectorDatum = reflector[elementID]
      if elementID == 0 {
        reflectorDatum = 1
      } else {
        reflectorDatum /= oldSubdiagonal - newSubdiagonal
      }
      reflectorDatum *= tau.squareRoot()
      
      if nanPresent {
        reflectorDatum = 0
      }
      reflector[elementID] = reflectorDatum
    }
    
    return reflector
  }
  
  // Applies the reflector to the trailing submatrix.
  private mutating func applyBulgeChase(
    reflector: [Float],
    startElementID: Int,
    endElementID: Int
  ) {
    let rangeCount = endElementID - startElementID
    
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
  }
}
