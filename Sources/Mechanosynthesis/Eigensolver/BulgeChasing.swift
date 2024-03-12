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
    let range = startElementID..<endElementID
    let indexOffset = startElementID
    let rangeCount = endElementID - startElementID
    
    // Load the row into the cache.
    var reflector = [Float](repeating: 0, count: rangeCount)
    for elementID in range {
      let address = vectorID * problemSize + elementID
      let matrixDatum = matrix[address]
      reflector[elementID - indexOffset] = matrixDatum
    }
    
    // Take the norm of the vector.
    var norm: Float = .zero
    for elementID in range {
      let reflectorDatum = reflector[elementID - indexOffset]
      norm += reflectorDatum * reflectorDatum
    }
    norm.formSquareRoot()
    
    // Predict the normalization factor.
    let oldSubdiagonal = reflector[startElementID - indexOffset]
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
    for elementID in range {
      var reflectorDatum = reflector[elementID - indexOffset]
      if elementID == startElementID {
        reflectorDatum = 1
      } else {
        reflectorDatum /= oldSubdiagonal - newSubdiagonal
      }
      reflectorDatum *= tau.squareRoot()
      
      if nanPresent {
        reflectorDatum = 0
      }
      reflector[elementID - indexOffset] = reflectorDatum
    }
    
    // Apply the reflector to the matrix, from both sides.
    for directionID in 0..<2 {
      // TODO: Simplify some inner loops by eliminatinating the subtraction of
      // 'indexOffset' whenever possible.
      // TODO: Predict which vectors will be affected by this reflector.
      for vectorID in 0..<problemSize {
        var vector = [Float](repeating: 0, count: rangeCount)
        
        if directionID == 0 {
          // Load the row into the cache.
          for elementID in range {
            let address = vectorID * problemSize + elementID
            let matrixDatum = matrix[address]
            vector[elementID - indexOffset] = matrixDatum
          }
        } else {
          // Load the column into the cache.
          for elementID in range {
            let address = elementID * problemSize + vectorID
            let matrixDatum = matrix[address]
            vector[elementID - indexOffset] = matrixDatum
          }
        }
        
        // Apply the reflector.
        var dotProduct: Float = .zero
        for elementID in range {
          let reflectorDatum = reflector[elementID - indexOffset]
          let vectorDatum = vector[elementID - indexOffset]
          dotProduct += reflectorDatum * vectorDatum
        }
        for elementID in range {
          let reflectorDatum = reflector[elementID - indexOffset]
          var vectorDatum = vector[elementID - indexOffset]
          vectorDatum -= reflectorDatum * dotProduct
          vector[elementID - indexOffset] = vectorDatum
        }
        
        if directionID == 0 {
          // Store the row to main memory.
          for elementID in range {
            let address = vectorID * problemSize + elementID
            let vectorDatum = vector[elementID - indexOffset]
            matrix[address] = vectorDatum
          }
        } else {
          // Store the column to main memory.
          for elementID in range {
            let address = elementID * problemSize + vectorID
            let vectorDatum = vector[elementID - indexOffset]
            matrix[address] = vectorDatum
          }
        }
      }
    }
    
    // Store the reflector to main memory.
    return reflector
  }
}
