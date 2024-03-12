//
//  BulgeChasing.swift
//
//
//  Created by Philip Turner on 3/11/24.
//

struct BulgeReflector {
  var data: [Float]
  var indices: Range<Int>
}

extension Diagonalization {
  // Returns a sequence of reflectors.
  mutating func chaseBulges() -> [BulgeReflector] {
    // If the matrix is already in tridiagonal form, there is no work to do.
    guard blockSize > 1 else {
      return []
    }
    
    var bulgeReflectors: [BulgeReflector] = []
    
    for sweepID in 0..<max(0, problemSize - 2) {
      let startVectorID = sweepID + 1
      var endVectorID = sweepID + blockSize + 1
      endVectorID = min(endVectorID, problemSize)
      guard endVectorID - startVectorID > 1 else {
        fatalError("Generated empty Householder transform.")
      }
      
      let data = applyBulgeChase(
        sweepID: sweepID, vectorID: sweepID,
        startElementID: startVectorID, endElementID: endVectorID)
      let range = startVectorID..<endVectorID
      let reflector = BulgeReflector(
        data: data, indices: range)
      bulgeReflectors.append(reflector)
      
      var operationID = 1
      while true {
        let offset = operationID * blockSize
        let startColumnID = (sweepID - blockSize + 1) + offset
        let startVectorID = (sweepID + 1) + offset
        var endVectorID = (sweepID + blockSize + 1) + offset
        endVectorID = min(endVectorID, problemSize)
        
        if endVectorID - startVectorID > 1 {
          let data = applyBulgeChase(
            sweepID: sweepID, vectorID: startColumnID,
            startElementID: startVectorID, endElementID: endVectorID)
          let range = startVectorID..<endVectorID
          let reflector = BulgeReflector(
            data: data, indices: range)
          bulgeReflectors.append(reflector)
        } else {
          break
        }
        operationID += 1
      }
    }
    
    return bulgeReflectors
  }
  
  // Returns a Householder reflector as an array allocation.
  //
  // Applies the Householder reflector to the entire matrix. This isn't very
  // efficient, but it's correct, which we want for initial debugging.
  private mutating func applyBulgeChase(
    sweepID: Int,
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
    
    // Allocate cache memory for the recipient of the Householder reflector.
    var vectorCache = [Float](repeating: 0, count: rangeCount)
    
    // Apply the reflector to the matrix, from both sides.
    for directionID in 0..<2 {
      let startApplicationID: Int = max(startElementID - blockSize, 0)
      let endApplicationID: Int = min(endElementID + blockSize, problemSize)
      for vectorID in startApplicationID..<endApplicationID {
        if directionID == 0 {
          // Load the row into the cache.
          for cacheElementID in 0..<rangeCount {
            let memoryElementID = startElementID + cacheElementID
            let matrixAddress = vectorID * problemSize + memoryElementID
            let matrixDatum = matrix[matrixAddress]
            vectorCache[cacheElementID] = matrixDatum
          }
        } else {
          // Load the column into the cache.
          for cacheElementID in 0..<rangeCount {
            let memoryElementID = startElementID + cacheElementID
            let matrixAddress = memoryElementID * problemSize + vectorID
            let matrixDatum = matrix[matrixAddress]
            vectorCache[cacheElementID] = matrixDatum
          }
        }
        
        // Apply the reflector.
        var dotProduct: Float = .zero
        for cacheElementID in 0..<rangeCount {
          let reflectorDatum = reflector[cacheElementID]
          let vectorDatum = vectorCache[cacheElementID]
          dotProduct += reflectorDatum * vectorDatum
        }
        for cacheElementID in 0..<rangeCount {
          let reflectorDatum = reflector[cacheElementID]
          var vectorDatum = vectorCache[cacheElementID]
          vectorDatum -= reflectorDatum * dotProduct
          vectorCache[cacheElementID] = vectorDatum
        }
        
        if directionID == 0 {
          // Store the row to main memory.
          for cacheElementID in 0..<rangeCount {
            let memoryElementID = startElementID + cacheElementID
            let matrixAddress = vectorID * problemSize + memoryElementID
            let vectorDatum = vectorCache[cacheElementID]
            matrix[matrixAddress] = vectorDatum
          }
        } else {
          // Store the column to main memory.
          for cacheElementID in 0..<rangeCount {
            let memoryElementID = startElementID + cacheElementID
            let matrixAddress = memoryElementID * problemSize + vectorID
            let vectorDatum = vectorCache[cacheElementID]
            matrix[matrixAddress] = vectorDatum
          }
        }
      }
    }
    
    // Store the reflector to main memory.
    return reflector
  }
}
