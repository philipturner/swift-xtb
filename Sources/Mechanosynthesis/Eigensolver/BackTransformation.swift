//
//  BackTransformation.swift
//
//
//  Created by Philip Turner on 3/11/24.
//

extension Diagonalization {
  mutating func backTransform(
    bulgeChasingReflectors: [BulgeReflector]
  ) {
    // Back-transform the eigenvectors.
    for vectorID in 0..<problemSize {
      // Load the vector into the cache.
      var vector = [Float](repeating: 0, count: problemSize)
      for elementID in 0..<problemSize {
        let address = vectorID * problemSize + elementID
        vector[elementID] = eigenvectors[address]
      }
      
      for reflectorID in bulgeChasingReflectors.indices.reversed() {
        // Load the reflector into the cache.
        let reflector = bulgeChasingReflectors[reflectorID]
        let indexOffset = reflector.indices.lowerBound
        
        // Apply the reflector.
        var dotProduct: Float = .zero
        for dataElementID in reflector.data.indices {
          let reflectorDatum = reflector.data[dataElementID]
          let vectorDatum = vector[indexOffset + dataElementID]
          dotProduct += reflectorDatum * vectorDatum
        }
        for dataElementID in reflector.data.indices {
          let reflectorDatum = reflector.data[dataElementID]
          var vectorDatum = vector[indexOffset + dataElementID]
          vectorDatum -= reflectorDatum * dotProduct
          vector[indexOffset + dataElementID] = vectorDatum
        }
      }
      
      // Store the vector to main memory.
      for elementID in 0..<problemSize {
        let address = vectorID * problemSize + elementID
        eigenvectors[address] = vector[elementID]
      }
    }
  }
  
  mutating func backTransform(
    bandFormReflectors: [Float]
  ) {
    // Back-transform the eigenvectors.
    for vectorID in 0..<problemSize {
      // Load the vector into the cache.
      var vector = [Float](repeating: 0, count: problemSize)
      for elementID in 0..<problemSize {
        let address = vectorID * problemSize + elementID
        vector[elementID] = eigenvectors[address]
      }
      
      
      for reflectorID in (0..<problemSize).reversed() {
        // Load the reflector into the cache.
        var reflector = [Float](repeating: 0, count: problemSize)
        for elementID in 0..<problemSize {
          let address = reflectorID * problemSize + elementID
          reflector[elementID] = bandFormReflectors[address]
        }
        
        // Apply the reflector.
        var dotProduct: Float = .zero
        for elementID in 0..<problemSize {
          dotProduct += reflector[elementID] * vector[elementID]
        }
        for elementID in 0..<problemSize {
          vector[elementID] -= reflector[elementID] * dotProduct
        }
      }
      
      // Store the vector to main memory.
      for elementID in 0..<problemSize {
        let address = vectorID * problemSize + elementID
        eigenvectors[address] = vector[elementID]
      }
    }
  }
}
