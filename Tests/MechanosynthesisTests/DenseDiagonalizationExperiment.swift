import XCTest
import Numerics

// All matrices in this experiment are assumed to be row-major as well. That
// may change when utilizing the AMX, because NT multiplications can achieve
// higher performance. For the final optimization, GPU acceleration seems more
// cost-effective than exploiting symmetry.
//
// The individual wavefunctions are rows of the eigenvector matrix.
final class DenseDiagonalizationExperiment: XCTestCase {
  // MARK: - Algorithms
  
  // Orthonormalizes the matrix in mixed precision, as a reference
  // implementation for more approximate methods.
  static func gramSchmidtOrthonormalize(
    matrix: inout [Float], n: Int
  ) {
    func normalize(electronID: Int) {
      var norm: Double = .zero
      for cellID in 0..<n {
        let value = matrix[electronID * n + cellID]
        norm += Double(value * value)
      }
      
      let normalizationFactor = 1 / norm.squareRoot()
      for cellID in 0..<n {
        var value = Double(matrix[electronID * n + cellID])
        value *= normalizationFactor
        matrix[electronID * n + cellID] = Float(value)
      }
    }
    
    for electronID in 0..<n {
      // Normalize the vectors before taking dot products.
      normalize(electronID: electronID)
    }
    
    for electronID in 0..<n {
      // Determine the magnitude of components parallel to previous vectors.
      var dotProducts = [Double](repeating: .zero, count: n)
      for neighborID in 0..<electronID {
        var dotProduct: Double = .zero
        for cellID in 0..<n {
          let value1 = matrix[electronID * n + cellID]
          let value2 = matrix[neighborID * n + cellID]
          dotProduct += Double(value1) * Double(value2)
        }
        dotProducts[neighborID] = dotProduct
      }
      
      // Subtract all components parallel to previous vectors.
      for cellID in 0..<n {
        var value1 = Double(matrix[electronID * n + cellID])
        for neighborID in 0..<electronID {
          let dotProduct = dotProducts[neighborID]
          let value2 = matrix[neighborID * n + cellID]
          value1 -= dotProduct * Double(value2)
        }
        matrix[electronID * n + cellID] = Float(value1)
      }
      
      // Rescale the orthogonal component to unit vector length.
      normalize(electronID: electronID)
    }
  }
  
  // This currently runs on single-core CPU, without AMX acceleration.
  // - returns: Conservative estimate of the maximum residual.
  //
  // Try modifying the orthogonalization algorithm to accelerate it:
  //
  // Perform something like panel factorization, where the matrix is divide into
  // an octree. Each sequential step unlocks a larger panel that also uses more
  // parallelism. For the small sequential parts along the diagonal, something
  // akin to force-based orthogonalization could work.
  //
  // What about a hybrid algorithm that annihilates only the greatest sources
  // of error with sequential Gram-Schmidt?
  static func parallelOrthogonalizeIteration(
    matrix: inout [Float], n: Int
  ) -> Float {
    // Generate the weight matrix (dot product matrix).
    var W = [Float](repeating: .zero, count: n * n)
    for electronID in 0..<n {
      for neighborID in 0..<n {
        var dotProduct: Float = .zero
        for cellID in 0..<n {
          let value1 = matrix[electronID * n + cellID]
          let value2 = matrix[neighborID * n + cellID]
          dotProduct += value1 * value2
        }
        W[electronID * n + neighborID] = dotProduct
      }
    }
    var maxW: Float = .zero
    for electronID in 0..<n {
      for neighborID in 0..<electronID {
        let value = W[electronID * n + neighborID]
        maxW = max(maxW, value.magnitude)
      }
    }
    
    // Generate the force matrix (approximation to the action of Gram-Schmidt).
    /*
     W <- lower triangle (W) O(n^2)
     F = -W * Ψ              O(n^3)
     */
    var F = [Float](repeating: .zero, count: n * n)
    for electronID in 0..<n {
      for neighborID in 0..<n {
        var weight = W[electronID * n + neighborID]
        if neighborID < electronID {
          // pass
        } else {
          weight = 0
        }
        for cellID in 0..<n {
          var forceValue = F[electronID * n + cellID]
          var neighborValue = matrix[neighborID * n + cellID]
          forceValue -= weight * neighborValue
          F[electronID * n + cellID] = forceValue
        }
      }
    }
    
    // Construct the timestep values.
    /*
     l = vector-wise ||F||   O(n^2)
     construct λ (l)         O(n^2)
     */
    var l = [Float](repeating: .zero, count: n)
    for electronID in 0..<n {
      var forceLength: Float = .zero
      for cellID in 0..<n {
        let value = F[electronID * n + cellID]
        forceLength += value * value
      }
      forceLength.formSquareRoot()
      l[electronID] = forceLength
    }
    var λ = [Float](repeating: .zero, count: n * n)
    for electronID in 0..<n {
      for neighborID in 0..<n {
        let length1 = λ[electronID]
        let length2 = λ[neighborID]
        let maxLength = max(length1, length2)
        
        var timeStep: Float = 1
        if maxLength > timeStep {
          timeStep /= maxLength
        }
        λ[electronID * n + neighborID] = timeStep
      }
    }
    
    // Apply the damped force.
    /*
     F = -λW * Ψ             O(n^3)
     Ψ <- Ψ + F              O(n^2)
     */
    F = [Float](repeating: .zero, count: n * n)
    for electronID in 0..<n {
      for neighborID in 0..<n {
        var weight = W[electronID * n + neighborID]
        if neighborID < electronID {
          // pass
        } else {
          weight = 0
        }
        weight *= λ[electronID * n + neighborID]
        for cellID in 0..<n {
          var forceValue = F[electronID * n + cellID]
          var neighborValue = matrix[neighborID * n + cellID]
          forceValue -= weight * neighborValue
          F[electronID * n + cellID] = forceValue
        }
      }
    }
    
    for entryID in F.indices {
      matrix[entryID] += F[entryID]
    }
    
    return maxW
  }
  
  // MARK: - Tests
  
  func testGramSchmidt() throws {
    // Generate a pre-determined matrix, test that the current version of
    // Gram-Schmidt normalizes them.
    var matrix: [Float] = [
      7, 6, 5, 4, 3, 2, 1,
      6, 7, 5, 4, 3, 2, 1,
      2, -1, -1, 1, 2, 6, 8,
      0.1, 0.2, 0.5, 0.5, 0.1, 0.2, 0.5,
      -0.1, 0.2, -0.5, 0.5, -0.1, 0.2, -0.5,
      -1, -2, -3, -5, -7, -9, -10,
      69, 23, 9, -48, 7, 1, 9,
    ]
    guard matrix.count == 49 else {
      fatalError("Not a 7x7 matrix.")
    }
    Self.gramSchmidtOrthonormalize(matrix: &matrix, n: 7)
    
    var overlapMatrix = [Float](repeating: .zero, count: 49)
    for electronID in 0..<7 {
      for neighborID in 0..<7 {
        var overlap: Float = .zero
        for cellID in 0..<7 {
          let value1 = matrix[electronID * 7 + cellID]
          let value2 = matrix[neighborID * 7 + cellID]
          overlap += value1 * value2
        }
        overlapMatrix[electronID * 7 + neighborID] = overlap
      }
    }
    
    for electronID in 0..<7 {
      for neighborID in 0..<7 {
        let value = overlapMatrix[electronID * 7 + neighborID]
        if electronID == neighborID {
          XCTAssertEqual(value, 1.0, accuracy: 1e-5)
        } else {
          XCTAssertEqual(value, 0.0, accuracy: 1e-5)
        }
      }
    }
    
    print()
    for electronID in 0..<7 {
      for neighborID in 0..<7 {
        let value = matrix[electronID * 7 + neighborID]
//        if electronID == neighborID {
//          XCTAssertEqual(value, 1.0, accuracy: 1e-5)
//        } else {
//          XCTAssertEqual(value, 0.0, accuracy: 1e-5)
//        }
        print(value, terminator: ", ")
      }
      print()
    }
  }
  
  // A self-contained test for the parallel orthogonalize algorithm. This is
  // almost identical to the Gram-Schmidt test.
  func testParallelOrthogonalize() throws {
    // Generate a pre-determined matrix, test that the current version of
    // Gram-Schmidt normalizes them.
    var matrix: [Float] = [
      7, 6, 5, 4, 3, 2, 1,
      6, 7, 5, 4, 3, 2, 1,
      2, -1, -1, 1, 2, 6, 8,
      0.1, 0.2, 0.5, 0.5, 0.1, 0.2, 0.5,
      -0.1, 0.2, -0.5, 0.5, -0.1, 0.2, -0.5,
      -1, -2, -3, -5, -7, -9, -10,
      69, 23, 9, -48, 7, 1, 9,
    ]
    guard matrix.count == 49 else {
      fatalError("Not a 7x7 matrix.")
    }
    func normalize(electronID: Int) {
      var norm: Double = .zero
      for cellID in 0..<7 {
        let value = matrix[electronID * 7 + cellID]
        norm += Double(value * value)
      }
      
      let normalizationFactor = 1 / norm.squareRoot()
      for cellID in 0..<7 {
        var value = Double(matrix[electronID * 7 + cellID])
        value *= normalizationFactor
        matrix[electronID * 7 + cellID] = Float(value)
      }
    }
    
    for electronID in 0..<7 {
      normalize(electronID: electronID)
    }
    
    print()
    for _ in 0..<10 {
      let residual = Self.parallelOrthogonalizeIteration(matrix: &matrix, n: 7)
      for electronID in 0..<7 {
        normalize(electronID: electronID)
      }
      print("residual:", residual)
      print("overlap:")
      
      var overlapMatrix = [Float](repeating: .zero, count: 49)
      for electronID in 0..<7 {
        for neighborID in 0..<7 {
          var overlap: Float = .zero
          for cellID in 0..<7 {
            let value1 = matrix[electronID * 7 + cellID]
            let value2 = matrix[neighborID * 7 + cellID]
            overlap += value1 * value2
          }
          overlapMatrix[electronID * 7 + neighborID] = overlap
        }
      }
      
      for electronID in 0..<7 {
        print("-", terminator: " ")
        for neighborID in 0..<7 {
          let value = overlapMatrix[electronID * 7 + neighborID]
          print(String(format: "%.3f", value.magnitude), terminator: ", ")
        }
        print()
      }
    }
    
    var overlapMatrix = [Float](repeating: .zero, count: 49)
    for electronID in 0..<7 {
      for neighborID in 0..<7 {
        var overlap: Float = .zero
        for cellID in 0..<7 {
          let value1 = matrix[electronID * 7 + cellID]
          let value2 = matrix[neighborID * 7 + cellID]
          overlap += value1 * value2
        }
        overlapMatrix[electronID * 7 + neighborID] = overlap
      }
    }
    
    print()
    for electronID in 0..<7 {
      for neighborID in 0..<7 {
        let value = matrix[electronID * 7 + neighborID]
//        if electronID == neighborID {
//          XCTAssertEqual(value, 1.0, accuracy: 1e-5)
//        } else {
//          XCTAssertEqual(value, 0.0, accuracy: 1e-5)
//        }
        print(value, terminator: ", ")
      }
      print()
    }
  }
}
