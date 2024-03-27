//
//  BulgeChasing.swift
//
//
//  Created by Philip Turner on 3/11/24.
//

import Accelerate

struct BulgeSweep {
  var data: [Float]
}

extension Diagonalization {
  // TODO: Store the band reflectors in a compact matrix. Change 'indices' to
  // describe a range in the matrix. Then, predict the values of 'indices' on
  // the fly, instead of storing them explicitly.
  //
  // The very first incremental change can simply be aggregating the
  // transforms within a specific sweep. After that, it is a minor change to
  // redirect the start of the vector to a compact matrix.
  //
  // This optimization does not include batching of Householder transforms.
  // However, the internal storage format could be optimized to improve the
  // performance of batching. Store transforms with the first element at
  // column[0], instead of along the diagonal.
  
  // Returns a sequence of reflectors.
  mutating func chaseBulges() -> [BulgeSweep] {
    // If the matrix is already in tridiagonal form, there is no work to do.
    guard blockSize > 1 else {
      return []
    }
    
    var bulgeSweeps: [BulgeSweep] = []
    
    for sweepID in 0..<max(0, problemSize - 2) {
      let startVectorID = sweepID + 1
      var endVectorID = sweepID + blockSize + 1
      endVectorID = min(endVectorID, problemSize)
      guard endVectorID - startVectorID > 1 else {
        fatalError("Generated empty Householder transform.")
      }
      var sweepData = [Float](repeating: 0, count: problemSize)
      
      // Apply the very first reflector.
      let data = applyBulgeChase(
        sweepID: sweepID, vectorID: sweepID,
        startElementID: startVectorID, endElementID: endVectorID)
      
      let maxReflectorElementID = problemSize - sweepID - 1
      let startReflectorElementID = 0
      let endReflectorElementID = min(blockSize, maxReflectorElementID)
      let reflectorRange = startReflectorElementID..<endReflectorElementID
      for reflectorElementID in reflectorRange {
        let value = data[reflectorElementID - startReflectorElementID]
        sweepData[reflectorElementID] = value
      }
      
      // Apply the remaining reflectors.
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
          
          let maxReflectorElementID = problemSize - sweepID - 1
          let startReflectorElementID = operationID * blockSize
          let endReflectorElementID = min(
            startReflectorElementID + blockSize, maxReflectorElementID)
          let reflectorRange = startReflectorElementID..<endReflectorElementID
          for reflectorElementID in reflectorRange {
            let value = data[reflectorElementID - startReflectorElementID]
            sweepData[reflectorElementID] = value
          }
        } else {
          break
        }
        operationID += 1
      }
      
      let sweep = BulgeSweep(data: sweepData)
      bulgeSweeps.append(sweep)
    }
    
    return bulgeSweeps
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
    
    // Apply the reflector to the matrix, from the left.
    matrix.withContiguousMutableStorageIfAvailable {
      let startApplicationID: Int = max(startElementID - blockSize, 0)
      let endApplicationID: Int = min(endElementID + blockSize, problemSize)
      let dotProductCount = endApplicationID - startApplicationID
      var dotProducts = [Float](repeating: 0, count: dotProductCount)
      let matrixBaseAddress: Int = startApplicationID * problemSize + startElementID
      let matrix = $0.baseAddress! + matrixBaseAddress
      
      #if false
      for m in 0..<dotProductCount {
        var dotProduct: Float = .zero
        for n in 0..<1 {
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
          var dotProduct: Float = .zero
          for k in 0..<1 {
            let lhsValue = dotProducts[m * 1 + k]
            let rhsValue = reflector[n * 1 + k]
            dotProduct += lhsValue * rhsValue
          }
          matrix[m * problemSize + n] -= dotProduct
        }
      }
      #elseif false
      do {
        var TRANSA = CChar(Character("N").asciiValue!)
        var TRANSB = CChar(Character("T").asciiValue!)
        var M = Int32(rangeCount)
        var N = Int32(dotProductCount)
        var K = Int32(1)
        var ALPHA = Float(-1)
        var LDA = Int32(rangeCount)
        var BETA = Float(1)
        var LDB = Int32(dotProductCount)
        var LDC = Int32(problemSize)
        sgemm_(
          &TRANSA, &TRANSB, &M, &N, &K, &ALPHA, reflector, &LDA,
          dotProducts, &LDB, &BETA, matrix, &LDC)
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
    do {
      let startApplicationID: Int = max(startElementID - blockSize, 0)
      let endApplicationID: Int = min(endElementID + blockSize, problemSize)
      for vectorID in startApplicationID..<endApplicationID {
        var dotProduct: Float = .zero
        for cacheElementID in 0..<rangeCount {
          let memoryElementID = startElementID + cacheElementID
          let matrixAddress = memoryElementID * problemSize + vectorID
          let matrixDatum = matrix[matrixAddress]
          
          let reflectorDatum = reflector[cacheElementID]
          dotProduct += reflectorDatum * matrixDatum
        }
        for cacheElementID in 0..<rangeCount {
          let memoryElementID = startElementID + cacheElementID
          let matrixAddress = memoryElementID * problemSize + vectorID
          var matrixDatum = matrix[matrixAddress]
          
          let reflectorDatum = reflector[cacheElementID]
          matrixDatum -= reflectorDatum * dotProduct
          matrix[matrixAddress] = matrixDatum
        }
      }
    }
    
    // Store the reflector to main memory.
    return reflector
  }
}
