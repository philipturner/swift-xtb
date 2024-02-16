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
  
  static func fastOrthogonalize(
    _ phi: [[Real]],
    priorities: [Real]? = nil,
    d3r: Real
  ) -> [[Real]] {
    var Ψ: [[Real]] = []
    for electronID in phi.indices {
      var currentΨ = phi[electronID]
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
          
          var weight: Real = 0.5
          if let priorities {
            let currentPriority = priorities[electronID]
            let otherPriority = priorities[otherID]
            weight = (otherPriority > currentPriority) ? 1 : 0
          }
          for cellID in force.indices {
            force[cellID] -= weight * dotProduct * otherΨ[cellID]
          }
        }
        
        let squareLength = integral(force, force, d3r: d3r)
        forceLengths.append(squareLength.squareRoot())
      }
      if maxDotProduct < 1e-7, iterationID > 15 {
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
          
          let currentLength = forceLengths[electronID].magnitude
          let otherLength = forceLengths[otherID].magnitude
          let maxLength = max(currentLength, otherLength)
          
          // The threshold of 0.5 seems to work for almost every random matrix.
          // I only saw one example where it failed to converge. It might be
          // caused by an especially ill-conditioned data set, which can be
          // avoided through other means.
          var weight: Real = 0.5
          if let priorities {
            let currentPriority = priorities[electronID]
            let otherPriority = priorities[otherID]
            weight = (otherPriority > currentPriority) ? 1 : 0
          }
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
  
  // Goal:
  // - Begin with the conventional Gram-Schmidt orthogonalizer. Then, switch to
  //   the iterative one and check that results are correct.
  // - For this test, use the eigenvalue as a "mass" for the "orthogonalization
  //   force". See how that affects stability and convergence rate.
  //
  // While running this experiment, I found an issue with the diagonalization
  // failing to diverge. I used the same approach that worked for 3x3 matrices
  // in MM4RigidBody, but it wouldn't converge to a reasonable wavefunction. In
  // an attempt to fix this, I used the steepest descent eigensolver from INQ.
  //
  // It worked for the one-electron case, but was mildly slow. For the multi-
  // electron case, many electrons diverged with positive kinetic energy. To fix
  // this, I made an extremely large, flat potential. Convergence was glacial
  // with Gram-Schmidt orthogonalization. The custom "fast" orthogonalizer
  // failed to converge at all. However, I don't think the fault is in the
  // orthogonalizer. It looks like something deeper.
  //
  // After debugging the particle in a finite potential well, and rewriting this
  // experiment from scratch, it worked. The "fast orthogonalizer" did converge
  // just as fast as the Gram-Schmidt orthogonalizer. A modification had to be
  // made though: bias the weights of the forces to either 0 or 1. So it has
  // the same priority order behavior as GS, just without the serial data
  // dependency between elements.
  
  func testOrthogonalizationExperiment() throws {
    // Hamiltonian: [-0.5 ∇^2 + (-1) / |x|] Ψ = E Ψ
    // Grid bounds: 0 Bohr - 3 Bohr
    // Grid spacing: 0.03 Bohr
    // 100 cells, 10 electrons
    //
    // Boundary conditions:
    // - wavefunction left of 0 Bohr equals wavefunction immediately right
    // - wavefunction right of 3 Bohr is zero
    
    let numCells = 100
    let numElectrons = 10
    let h: Real = 0.6 // the original 0.03 caused convergence issues
    let coreCharge: Real = 1
    
    func hamiltonian(_ Ψ: [[Real]]) -> [[Real]] {
      var output: [[Real]] = []
      for electronID in Ψ.indices {
        var outputVector: [Real] = []
        let currentΨ = Ψ[electronID]
        
        for cellID in currentΨ.indices {
          let filter: [Real] = [1, -2, 1]
          
          var value = filter[1] * currentΨ[cellID]
          if cellID < numCells - 1 {
            value += filter[2] * currentΨ[cellID + 1]
          }
          if cellID > 0 {
            value += filter[0] * currentΨ[cellID - 1]
          } else {
            value += filter[0] * currentΨ[cellID]
          }
          value *= -0.5
          value /= (h * h)
          
          let x = (Real(cellID) + 0.5) * h
          let v = -coreCharge / x
          value += v * currentΨ[cellID]
          
          outputVector.append(value)
        }
        output.append(outputVector)
      }
      return output
    }
    
    func dot(_ lhs: [Real], _ rhs: [Real]) -> Real {
      guard lhs.count == rhs.count else {
        fatalError("Unequal counts.")
      }
      var sum: Double = .zero
      for cellID in lhs.indices {
        sum += Double(lhs[cellID] * rhs[cellID] * h)
      }
      
      return Real(sum)
    }
    
    func normalize(_ x: [Real]) -> [Real] {
      let norm = dot(x, x)
      return x.map {
        $0 / norm.squareRoot()
      }
    }
    
    var Ψ: [[Real]] = []
    for electronID in 0..<numElectrons {
      var outputVector: [Real] = []
      let angularFrequency = Real(electronID) * 2 * Real.pi
      
      for cellID in 0..<numCells {
        let x = (Real(cellID) + 0.5) * h
        let value = exp(-x / 10) * cos(angularFrequency * x / 10)
        outputVector.append(value)
      }
      
      outputVector = normalize(outputVector)
      Ψ.append(outputVector)
    }
    
    func eigensolver(
      hamiltonian: ([[Real]]) -> [[Real]],
      phi originalΨ: [[Real]]
    ) -> [[Real]] {
      var Ψ = originalΨ
      var HΨ = hamiltonian(originalΨ)
      
      // Try just 1 eigensolver step to demonstrate a non-negligible
      // improvement to stability.
      for _ in 0..<5 {
        let E = zip(Ψ, HΨ).map(dot(_:_:))
        let norm = zip(Ψ, Ψ).map(dot(_:_:))
        
        var r: [[Real]] = []
        for electronID in Ψ.indices {
          let currentΨ = Ψ[electronID]
          let currentHΨ = HΨ[electronID]
          var outputVector: [Real] = []
          let evnorm = E[electronID] / norm[electronID]
          
          for cellID in currentΨ.indices {
            let value = currentHΨ[cellID] - evnorm * currentΨ[cellID]
            outputVector.append(value)
          }
          r.append(outputVector)
        }
        
        let Hr = hamiltonian(r)
        let rr = zip(r, r).map(dot(_:_:))
        let Ψr = zip(Ψ, r).map(dot(_:_:))
        let rHr = zip(r, Hr).map(dot(_:_:))
        let ΨHr = zip(Ψ, Hr).map(dot(_:_:))
        
        var λ: [Real] = []
        for electronID in 0..<Ψ.count {
          let m0 = rr[electronID]
          let m1 = Ψr[electronID]
          let m2 = rHr[electronID]
          let m3 = ΨHr[electronID]
          let m4 = E[electronID]
          let m5 = norm[electronID]
          
          let ca = (m0 * m3 - m2 * m1)
          let cb = (m5 * m2 - m4 * m0)
          let cc = (m4 * m1 - m3 * m5)
          let determinant = cb + (cb * cb - 4 * ca * cc).squareRoot()
          if determinant.magnitude < 1e-15 {
            // This happens if we are perfectly converged.
            λ.append(0)
          } else {
            λ.append(2 * cc / determinant)
          }
        }
        
        for electronID in Ψ.indices {
          for cellID in Ψ[electronID].indices {
            let lambda = λ[electronID]
            Ψ[electronID][cellID] += lambda * r[electronID][cellID]
            HΨ[electronID][cellID] += lambda * Hr[electronID][cellID]
          }
        }
      }
      
      for electronID in Ψ.indices {
        Ψ[electronID] = normalize(Ψ[electronID])
      }
      
      // We proved that "fast orthogonalization" works just as well as regular
      // orthogonalization. However, we'll use the original orthogonalization
      // function for this unit test. Fast orthogonalization is not optimized
      // in its current form, and this test converges very slowly.
      Ψ = OrthogonalizeTests.orthogonalize(Ψ, d3r: h)
//      let priorities = (0..<numElectrons).map { Real(-$0) }
//      Ψ = OrthogonalizeTests.fastOrthogonalize(
//        Ψ, priorities: priorities, d3r: h)
      return Ψ
    }
    
    let numIterations = 300
    for iterationID in 0...numIterations {
      /*
       solvers::steepest_descent(laplacian, identity, phi);
       
       auto residual = laplacian(phi);
       auto eigenvalues = operations::overlap_diagonal(phi, residual);
       operations::shift(eigenvalues, phi, residual, -1.0);
       auto normres = operations::overlap_diagonal(residual);
       
       for(int ivec = 0; ivec < phi.set_size(); ivec++){
       tfm::format(std::cout, " state %4d  evalue = %18.12f  res = %5.0e\n", ivec + 1, real(eigenvalues[ivec]), real(normres[ivec]));
       }
       */
      
      if OrthogonalizeTests.console {
        print("Iteration \(iterationID)")
      }
      
      if iterationID > 0 {
        Ψ = eigensolver(hamiltonian: hamiltonian(_:), phi: Ψ)
      }
      
      let HΨ = hamiltonian(Ψ)
      var E: [Real] = []
      var r: [[Real]] = []
      for electronID in 0..<numElectrons {
        let currentΨ = Ψ[electronID]
        let currentHΨ = HΨ[electronID]
        let currentE = dot(currentΨ, currentHΨ)
        E.append(currentE)
        
        var rVector: [Real] = []
        for cellID in 0..<numCells {
          let value = currentHΨ[cellID] - currentE * currentΨ[cellID]
          rVector.append(value)
        }
        r.append(rVector)
      }
      
      for electronID in 0..<numElectrons {
        let eigenvalue = E[electronID]
        let normres = dot(r[electronID], r[electronID])
        
        if iterationID == numIterations {
          XCTAssertLessThan(
            normres.magnitude, 5e-6,
            "Electron \(electronID) failed to converge.")
          
          let expectedEigenvalues: [Real] = [
            -2.4700130,
             -0.23662587,
             -0.08361732,
             -0.04222589,
             -0.025348555,
             -0.0155251175,
             -0.004280047,
             0.010540266,
             0.028675715,
             0.049906835,
          ]
          XCTAssertEqual(
            eigenvalue, expectedEigenvalues[electronID], accuracy: 1e-4,
            "Electron \(electronID) had the wrong eigenvalue.")
        }
        
        if OrthogonalizeTests.console {
          var output = ""
          output += "  state \(electronID)"
          output += "  evalue = \(eigenvalue)"
          output += "  res = \(normres)"
          print(output)
        }
      }
      
      if OrthogonalizeTests.console, iterationID == numIterations {
        for cellID in 0..<numCells {
          for electronID in 0..<numElectrons {
            let value = Ψ[electronID][cellID]
            var repr = String(format: "%.4f", value)
            if !repr.starts(with: "-") {
              repr = " " + repr
            }
            print("  \(repr),", separator: "", terminator: "")
          }
          print()
        }
      }
    }
  }
}
