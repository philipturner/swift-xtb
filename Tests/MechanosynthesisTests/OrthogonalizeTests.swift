import XCTest
import Mechanosynthesis
import Numerics

// Test an alternative to Cholesky decomposition that's more parallel.
final class OrthogonalizeTests: XCTestCase {
  typealias Real = Float
  
  static let console: Bool = false
  
  static func integral(_ lhs: [Real], _ rhs: [Real], d3r: Real) -> Real {
    precondition(lhs.count == rhs.count, "Vectors had different count.")
    var sum: Double = .zero
    for cellID in lhs.indices {
      sum += Double(lhs[cellID] * rhs[cellID] * d3r)
    }
    return Real(sum)
  }
  
  // Accepts a (not necessarily normalized) set of vectors and returns an
  // orthonormal set. Gives each vector the same weight.
  static func orthogonalize(_ phi: [[Real]], d3r: Real) -> [[Real]] {
    var orthogonalizedΨ: [[Real]] = []
    for electronID in phi.indices {
      var currentΨ = phi[electronID]
      var norm = integral(currentΨ, currentΨ, d3r: d3r)
      for cellID in currentΨ.indices {
        currentΨ[cellID] /= norm.squareRoot()
      }
      
      var dotProducts: [Real] = []
      for otherID in orthogonalizedΨ.indices {
        let otherΨ = orthogonalizedΨ[otherID]
        dotProducts.append(integral(currentΨ, otherΨ, d3r: d3r))
      }
      for otherID in orthogonalizedΨ.indices {
        let dotProduct = dotProducts[otherID]
        let otherΨ = orthogonalizedΨ[otherID]
        for cellID in currentΨ.indices {
          currentΨ[cellID] -= dotProduct * otherΨ[cellID]
        }
      }
      
      norm = integral(currentΨ, currentΨ, d3r: d3r)
      for cellID in currentΨ.indices {
        currentΨ[cellID] /= norm.squareRoot()
      }
      orthogonalizedΨ.append(currentΨ)
    }
    return orthogonalizedΨ
  }
  
  static func fastOrthogonalize(_ unnormalizedPhi: [[Real]], d3r: Real) -> [[Real]] {
    var Ψ: [[Real]] = []
    for electronID in unnormalizedPhi.indices {
      var currentΨ = unnormalizedPhi[electronID]
      let norm = integral(currentΨ, currentΨ, d3r: d3r)
      for cellID in currentΨ.indices {
        currentΨ[cellID] /= norm.squareRoot()
      }
      Ψ.append(currentΨ)
    }
    
    let maxIterations = 20
    for iterationID in 1...maxIterations {
      // Find how much each vector wants to change.
      var maxDotProduct: Real = .zero
      var forceLengths: [Real] = []
      for electronID in Ψ.indices {
        let currentΨ = Ψ[electronID]
        var force = Array(repeating: Real(0), count: currentΨ.count)
        for otherID in Ψ.indices where electronID != otherID {
          let otherΨ = Ψ[otherID]
          let dotProduct = integral(currentΨ, otherΨ, d3r: d3r)
          maxDotProduct = max(maxDotProduct, dotProduct.magnitude)
          
          let weight: Real = 0.5
          for cellID in force.indices {
            force[cellID] -= weight * dotProduct * otherΨ[cellID]
          }
        }
        
        let squareLength = integral(force, force, d3r: d3r)
        forceLengths.append(squareLength.squareRoot())
      }
      if maxDotProduct < 1e-7 {
        if console {
          print("Converged on iteration \(iterationID)")
        }
        break
      } else if iterationID == maxIterations {
        print("Failed to converge after \(maxIterations) iterations. Maximum dot product was \(maxDotProduct).")
        break
      } else {
        if console {
          print("iteration \(iterationID): \(maxDotProduct)")
        }
      }
      
      // Change the vectors.
      var orthogonalizedΨ: [[Real]] = []
      for electronID in Ψ.indices {
        var currentΨ = Ψ[electronID]
        var force = Array(repeating: Real(0), count: currentΨ.count)
        for otherID in Ψ.indices where electronID != otherID {
          let otherΨ = Ψ[otherID]
          let dotProduct = integral(currentΨ, otherΨ, d3r: d3r)
          
          let currentLength = forceLengths[electronID]
          let otherLength = forceLengths[otherID]
          let maxLength = max(currentLength, otherLength)
          
          // The threshold of 0.5 seems to work for almost every random matrix.
          // I only saw one example where it failed to converge. It might be
          // caused by an especially ill-conditioned data set, which can be
          // avoided through other means.
          var weight: Real = 0.5
          let threshold: Real = 0.5
          if maxLength > threshold {
            weight *= threshold / maxLength
          }
          for cellID in force.indices {
            force[cellID] -= weight * dotProduct * otherΨ[cellID]
          }
        }
        
        for cellID in force.indices {
          currentΨ[cellID] += force[cellID]
        }
        let norm = integral(currentΨ, currentΨ, d3r: d3r)
        for cellID in currentΨ.indices {
          currentΨ[cellID] /= norm.squareRoot()
        }
        orthogonalizedΨ.append(currentΨ)
      }
      
      // Prepare for the next iteration.
      Ψ = orthogonalizedΨ
    }
    
    return Ψ
  }
  
  // Outputs a symmetric matrix.
  static func overlap(_ phi: [[Real]], d3r: Real) -> [[Real]] {
    var output: [[Real]] = []
    for electronID in phi.indices {
      let currentΨ = phi[electronID]
      var outputVector: [Real] = []
      for otherID in phi.indices {
        let otherΨ = phi[otherID]
        let sum = integral(currentΨ, otherΨ, d3r: d3r)
        outputVector.append(sum)
      }
      output.append(outputVector)
    }
    return output
  }
  
  // Reproduce some of the tests from:
  // https://gitlab.com/npneq/inq/-/blob/master/src/operations/orthogonalize.hpp
  //
  // Use real numbers instead of complex numbers.
  func testDimension() throws {
    let spacing: Real = 0.44428829
    let gridWidth = Int((3.15 / spacing).rounded(.up))
    
    func createPhi(numPoints: Int) -> [[Real]] {
      var phi: [[Real]] = []
      for _ in 0..<numPoints {
        var phiVector: [Real] = []
        for _ in 0..<(gridWidth * gridWidth * gridWidth) {
          let value = Real.random(in: 0..<1)
          phiVector.append(value)
        }
        
        let norm = Self.integral(
          phiVector, phiVector, d3r: spacing * spacing * spacing)
        for cellID in 0..<(gridWidth * gridWidth * gridWidth) {
          phiVector[cellID] /= norm.squareRoot()
        }
        phi.append(phiVector)
      }
      return phi
    }
    
    func displayOverlapMatrix(phi: [[Real]]) {
      let overlap = Self.overlap(phi, d3r: spacing * spacing * spacing)
      
      if Self.console {
        print("overlap matrix:")
        for electronID in phi.indices {
          var output: String = ""
          for otherID in phi.indices {
            let value = overlap[electronID][otherID]
            let repr = "\(value)"
            if !repr.starts(with: "-") {
              output += " "
            }
            output += repr
            output += " "
          }
          print(output)
        }
      }
    }
    
    func checkOverlapMatrix(phi: [[Real]], accuracy: Real) {
      let overlap = Self.overlap(phi, d3r: spacing * spacing * spacing)
      for electronID in overlap.indices {
        var offDiagonalIsZero = true
        for otherID in overlap.indices {
          let value = overlap[electronID][otherID]
          if electronID == otherID {
            XCTAssertEqual(value, 1, accuracy: accuracy)
          } else {
            // Avoid printing an excessive number of errors to the console.
            if value.magnitude > accuracy {
              offDiagonalIsZero = false
            }
          }
        }
        XCTAssert(offDiagonalIsZero)
      }
    }
    
    func markNewTest() {
      if Self.console {
        print()
      }
    }
    
    // Dimension 3
    do {
      markNewTest()
      var phi = createPhi(numPoints: 3)
      displayOverlapMatrix(phi: phi)
      phi = Self.fastOrthogonalize(phi, d3r: spacing * spacing * spacing)
      displayOverlapMatrix(phi: phi)
      checkOverlapMatrix(phi: phi, accuracy: 1e-6)
    }
    
    // Dimension 7
    do {
      markNewTest()
      var phi = createPhi(numPoints: 7)
      displayOverlapMatrix(phi: phi)
      phi = Self.fastOrthogonalize(phi, d3r: spacing * spacing * spacing)
      displayOverlapMatrix(phi: phi)
      checkOverlapMatrix(phi: phi, accuracy: 1e-6)
    }
    
    // Dimension 7, Orthogonalize 15 Times
    do {
      markNewTest()
      var phi = createPhi(numPoints: 7)
      for _ in 0..<15 {
        phi = Self.fastOrthogonalize(phi, d3r: spacing * spacing * spacing)
      }
      checkOverlapMatrix(phi: phi, accuracy: 1e-6)
    }
    
    // Dimension 37
    do {
      markNewTest()
      var phi = createPhi(numPoints: 37)
      phi = Self.fastOrthogonalize(phi, d3r: spacing * spacing * spacing)
      checkOverlapMatrix(phi: phi, accuracy: 1e-6)
    }
  }
  
  // TODO: Run the orthogonalization experiment that was previously abandoned.
  // - Migrate the code from documentation to the test suite.
  // - For this test, use the eigenvalue as a "mass" for the "orthogonalization
  //   force". See how that affects stability and convergence rate.
  // - Begin with the conventional Gram-Schmidt orthogonalizer. Then, switch to
  //   the iterative one and check that results are correct.
}
