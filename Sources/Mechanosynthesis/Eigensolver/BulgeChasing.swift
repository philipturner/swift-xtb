//
//  BulgeChasing.swift
//
//
//  Created by Philip Turner on 3/11/24.
//

extension Diagonalization {
  // Returns a sequence of reflectors.
  mutating func chaseBulges() -> [[Float]] {
    // If the matrix is already in tridiagonal form, there is no work to do.
    guard blockSize > 1 else {
      return []
    }
    
    var bulgeReflectors: [[Float]] = []
    
    for sweepID in 0..<max(0, problemSize - 2) {
      let startVectorID = sweepID + 1
      var endVectorID = sweepID + blockSize + 1
      endVectorID = min(endVectorID, problemSize)
      guard endVectorID - startVectorID > 1 else {
        fatalError("Generated empty Householder transform.")
      }
      
      let reflector = applyBulgeChase(
        sweepID: sweepID, vectorID: sweepID,
        startElementID: startVectorID, endElementID: endVectorID)
      bulgeReflectors.append(reflector)
      
      var operationID = 1
      while true {
        let offset = operationID * blockSize
        let startColumnID = (sweepID - blockSize + 1) + offset
        let startVectorID = (sweepID + 1) + offset
        var endVectorID = (sweepID + blockSize + 1) + offset
        endVectorID = min(endVectorID, problemSize)
        
        if endVectorID - startVectorID > 1 {
          let reflector = applyBulgeChase(
            sweepID: sweepID, vectorID: startColumnID,
            startElementID: startVectorID, endElementID: endVectorID)
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
    // Load the row into the cache.
    var vector = [Float](repeating: 0, count: problemSize)
    for elementID in startElementID..<endElementID {
      let address = vectorID * problemSize + elementID
      vector[elementID] = matrix[address]
    }
    
    // Take the norm of the vector.
    var norm: Float = .zero
    for elementID in startElementID..<endElementID {
      norm += vector[elementID] * vector[elementID]
    }
    norm.formSquareRoot()
    
    // Modify the vector, turning it into a reflector.
    let oldSubdiagonal = vector[startElementID]
    let newSubdiagonal = norm * Float((oldSubdiagonal >= 0) ? -1 : 1)
    let tau = (newSubdiagonal - oldSubdiagonal) / newSubdiagonal
    for elementID in startElementID..<endElementID {
      var element = vector[elementID]
      if elementID == startElementID {
        element = 1
      } else {
        element /= oldSubdiagonal - newSubdiagonal
      }
      element *= tau.squareRoot()
      vector[elementID] = element
    }
    
    // Apply the reflector to the matrix, from both sides.
    let reflector = vector
    for directionID in 0..<2 {
      for vectorID in 0..<problemSize {
        var vector = [Float](repeating: 0, count: problemSize)
        if directionID == 0 {
          // Load the row into the cache.
          for elementID in 0..<problemSize {
            let address = vectorID * problemSize + elementID
            vector[elementID] = matrix[address]
          }
        } else {
          // Load the column into the cache.
          for elementID in 0..<problemSize {
            let address = elementID * problemSize + vectorID
            vector[elementID] = matrix[address]
          }
        }
        
        // Apply the reflector.
        var dotProduct: Float = .zero
        for elementID in 0..<problemSize {
          dotProduct += reflector[elementID] * vector[elementID]
        }
        for elementID in 0..<problemSize {
          vector[elementID] -= reflector[elementID] * dotProduct
        }
        
        if directionID == 0 {
          // Store the row to main memory.
          for elementID in 0..<problemSize {
            let address = vectorID * problemSize + elementID
            matrix[address] = vector[elementID]
          }
        } else {
          // Store the column to main memory.
          for elementID in 0..<problemSize {
            let address = elementID * problemSize + vectorID
            matrix[address] = vector[elementID]
          }
        }
      }
    }
    
    // Store the reflector to main memory.
    return reflector
  }
}
