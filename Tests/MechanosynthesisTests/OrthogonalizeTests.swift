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
  
  static func fastOrthogonalize(_ phi: [[Real]], d3r: Real) -> [[Real]] {
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
#if false
  func testOrthogonalizationExperiment() {
    // Hamiltonian: [-0.5 ∇^2 + (-1) / |x|] Ψ = E Ψ
    // Grid bounds: 0 Bohr - 3 Bohr
    // Grid spacing: 0.03 Bohr
    // 100 cells, 10 electrons
    //
    // Boundary conditions:
    // - wavefunction left of 0 Bohr equals wavefunction immediately right
    // - wavefunction right of 3 Bohr is zero
    
    var Ψ: [[Float]] = []
    let gridSpacing: Float = 0.03
    
    srand48(2021)
    for electronID in 0..<10 {
      var currentΨ: [Float] = []
      for cellID in 0..<100 {
        let x = gridSpacing * (Float(cellID) + 0.5)
        let cosineFrequency = Float(9 - electronID) / 2
        currentΨ.append(exp(-x) * cos(cosineFrequency * x))
//        currentΨ.append(cos(cosineFrequency * x))
      }
      
      var sum: Double = .zero
      for fragment in currentΨ {
        sum += Double(fragment * fragment * gridSpacing)
      }
      let normalizationFactor = Float((1 / sum).squareRoot())
      for cellID in 0..<100 {
        currentΨ[cellID] *= normalizationFactor
      }
      Ψ.append(currentΨ)
    }
    
    func formatString(_ value: Double) -> String {
      let repr = String(format: "%.3f", value)
      if !repr.starts(with: "-") {
        return " " + repr
      } else {
        return repr
      }
    }
    
    func reportState(_ Ψ: [[Float]]) {
      print()
      print(
        "| norm", "|",
        "  x   ", "|",
        "-∇^2/2", "|",
        "-1/|x|", "|",
        "Ψ H Ψ ")
      for electronID in 0..<10 {
        let currentΨ = Ψ[electronID]
        var observable1: Double = .zero
        var observableX: Double = .zero
        var observableKinetic: Double = .zero
        var observablePotential: Double = .zero
        var observableHamiltonian: Double = .zero
        
        for cellID in 0..<100 {
          let value = currentΨ[cellID]
          let x = gridSpacing * (Float(cellID) + 0.5)
          let leftValue = (cellID > 0) ? currentΨ[cellID - 1] : value
          let rightValue = (cellID < 99) ? currentΨ[cellID + 1] : 0
          
          let derivativeLeft = (value - leftValue) / gridSpacing
          let derivativeRight = (rightValue - value) / gridSpacing
          let laplacian = (derivativeRight - derivativeLeft) / gridSpacing
//          let corePotential = -0 / x.magnitude
          let corePotential = Float(-100)
          
          let d3r = gridSpacing
          observable1 += Double(value * 1 * value * d3r)
          observableX += Double(value * x * value * d3r)
          observableKinetic += Double(value * -0.5 * laplacian * d3r)
          observablePotential += Double(value * corePotential * value * d3r)
          
          let hamiltonian = -0.5 * laplacian + corePotential * value
          observableHamiltonian += Double(value * hamiltonian * d3r)
        }
        print(
          formatString(observable1), "|",
          formatString(observableX), "|",
          formatString(observableKinetic), "|",
          formatString(observablePotential), "|",
          formatString(observableHamiltonian))
      }
    }
    
    func reportWaveFunctions(_ Ψ: [[Float]]) {
      print()
      for cellID in 0..<100 {
        var output: String = ""
        for electronID in 0..<10 {
          let value = Ψ[electronID][cellID]
          output += formatString(Double(value)) + " "
        }
        print(output)
      }
    }
    
    func hamiltonian(_ Ψ: [Float]) -> [Float] {
      let currentΨ = Ψ
      var HΨ: [Float] = []
      for cellID in 0..<100 {
        let value = currentΨ[cellID]
        let x = gridSpacing * (Float(cellID) + 0.5)
        let leftValue = (cellID > 0) ? currentΨ[cellID - 1] : value
        let rightValue = (cellID < 99) ? currentΨ[cellID + 1] : 0
        
        let derivativeLeft = (value - leftValue) / gridSpacing
        let derivativeRight = (rightValue - value) / gridSpacing
        let laplacian = (derivativeRight - derivativeLeft) / gridSpacing
//        let corePotential = -0 / x.magnitude
        let corePotential = Float(-100)
        
        let hamiltonian = -0.5 * laplacian + corePotential * value
        HΨ.append(hamiltonian)
      }
      return HΨ
    }
    
    // Create another function for performing integrals, test it on the already
    // known values for energy.
    func integral(_ lhs: [Float], _ rhs: [Float]) -> Float {
      precondition(lhs.count == rhs.count, "Vectors had different count.")
      var sum: Double = .zero
      for cellID in lhs.indices {
        sum += Double(lhs[cellID] * rhs[cellID] * gridSpacing)
      }
      return Float(sum)
    }
    
    
    // This acts on one eigenvector at once. The orthogonalization must occur
    // outside of this function.
    func eigensolver(_ phi: [Float], maxResidual: inout Float) -> [Float] {
      var Ψ = phi
      var HΨ = hamiltonian(Ψ)
      
      for _ in 0..<5 {
        let E = integral(Ψ, HΨ)
        let norm = integral(Ψ, Ψ)
        let evnorm = E / norm
        
        var r: [Float] = []
        for cellID in 0..<100 {
          r.append(HΨ[cellID] - evnorm * Ψ[cellID])
        }
        
        let Hr = hamiltonian(r)
        let rr  = integral(r, r)
        let Ψr = integral(phi, r)
        let rHr = integral(r, Hr)
        let ΨHr = integral(phi, Hr)
        //      print(rr, Ψr, rHr, ΨHr)
        /*
         184.59229 -9.0687536e-08 71424.17 184.59229
         147.84608 3.2480602e-07 23955.428 147.84607
         224.95732 2.2965014e-07 95228.98 224.95726
         181.14989 4.800313e-07 36294.53 181.14996
         264.35968 1.4699713e-06 115933.73 264.35962
         210.91487 6.4571213e-07 47135.383 210.91513
         287.8893 1.69327e-07 119468.195 287.8893
         240.79008 -5.905458e-07 57412.242 240.79025
         303.56497 2.6445014e-07 114152.34 303.5647
         271.39746 1.9259278e-07 68358.49 271.39746
         */
        
        let m0 = rr
        let m1 = Ψr
        let m2 = rHr
        let m3 = ΨHr
        let m4 = E
        let m5 = norm
        
        let ca = m0 * m3 - m2 * m1
        let cb = m5 * m2 - m4 * m0
        let cc = m4 * m1 - m3 * m5
        //      print(ca, cb, cc)
        /*
         34074.32 72873.305 -184.5923
         21858.455 25227.379 -147.84607
         50605.758 97405.62 -224.95726
         32815.28 38169.742 -181.14996
         69885.85 118601.41 -264.35962
         44485.105 49124.12 -210.91513
         82880.234 121806.375 -287.88928
         57979.94 59005.496 -240.79025
         92151.58 115565.62 -303.5647
         73656.56 69041.266 -271.39746
         */
        
        let determinant = cb + (cb * cb - 4 * ca * cc).squareRoot()
        let λ = 2 * cc / determinant
        //      print(determinant, λ)
        /*
         145919.03 -0.0025300647
         50709.67 -0.0058310796
         195044.7 -0.0023067251
         76649.7 -0.004726697
         237513.95 -0.0022260556
         98628.766 -0.0042769494
         244003.9 -0.0023597104
         118482.32 -0.0040645767
         231614.34 -0.0026212945
         138659.2 -0.003914597
         */
//        print("  - λ = \(λ), rms_norm = \(rr.squareRoot())")
        maxResidual = max(maxResidual, rr.squareRoot())
        
        for cellID in 0..<100 {
          Ψ[cellID] += λ * r[cellID]
          HΨ[cellID] += λ * Hr[cellID]
        }
      }
      
      let norm = integral(Ψ, Ψ)
      for cellID in 0..<100 {
        Ψ[cellID] /= norm.squareRoot()
      }
      return Ψ
    }
    
    reportState(Ψ)
//    reportWaveFunctions(Ψ)
    
    for _ in 0..<30 {
      var maxResidual: Float = .zero
      for electronID in 0..<10 {
//        print("- electron \(electronID):")
        Ψ[electronID] = eigensolver(Ψ[electronID], maxResidual: &maxResidual)
      }
      Ψ = Self.orthogonalize(Ψ, d3r: gridSpacing)
      
//      reportState(Ψ)
      
      print("max residual: \(maxResidual)")
      
    }
    
    reportState(Ψ)
//    reportWaveFunctions(Ψ)
  }
  #endif
}
