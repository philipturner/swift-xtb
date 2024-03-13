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
    bandFormReflectors: [BandReflector]
  ) {
    for bandReflector in bandFormReflectors.reversed() {
      let panelReflectors = bandReflector.matrixV
      let T = bandReflector.matrixT
      
      // V^H A
      var VA = [Float](repeating: 0, count: problemSize * blockSize)
      for m in 0..<blockSize {
        for n in 0..<problemSize {
          var dotProduct: Float = .zero
          for k in 0..<problemSize {
            let lhsValue = panelReflectors[m * problemSize + k]
            let rhsValue = eigenvectors[n * problemSize + k]
            dotProduct += lhsValue * rhsValue
          }
          VA[n * blockSize + m] = dotProduct
        }
      }
      
      // T^H (V^H A)
      var TVA = [Float](repeating: 0, count: problemSize * blockSize)
      for m in 0..<blockSize {
        for n in 0..<problemSize {
          var dotProduct: Float = .zero
          for k in 0..<blockSize {
            let lhsValue = T[k * blockSize + m]
            let rhsValue = VA[n * blockSize + k]
            dotProduct += lhsValue * rhsValue
          }
          TVA[n * blockSize + m] = dotProduct
        }
      }
      
      // V (T^H V^H A)
      for m in 0..<problemSize {
        for n in 0..<problemSize {
          var dotProduct: Float = .zero
          for k in 0..<blockSize {
            let lhsValue = panelReflectors[k * problemSize + m]
            let rhsValue = TVA[n * blockSize + k]
            dotProduct += lhsValue * rhsValue
          }
          eigenvectors[n * problemSize + m] -= dotProduct
        }
      }
    }
  }
}
