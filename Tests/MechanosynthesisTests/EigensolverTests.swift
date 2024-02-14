import XCTest
import Mechanosynthesis
import Numerics

// Test the steepest descent algorithm from INQ.
final class EigensolverTests: XCTestCase {
  static let console: Bool = false
  
  // Reproduce the test from the shorter file:
  // https://gitlab.com/npneq/inq/-/blob/master/src/solvers/steepest_descent.hpp
  //
  // Compute everything in FP64.
  func testDiagonalMatrixDouble() throws {
    // Create row-major hamiltonian matrices.
    let numPoints = 100
    let numVectors = 12
    var identityMatrix: [[Complex<Double>]] = []
    var diagonalMatrix: [[Complex<Double>]] = []
    
    for i in 0..<numPoints {
      var identityRow: [Complex<Double>] = []
      var diagonalRow: [Complex<Double>] = []
      for j in 0..<numPoints {
        if i == j {
          identityRow.append(Complex(1))
          diagonalRow.append(Complex(Double(i) + 1))
        } else {
          identityRow.append(.zero)
          diagonalRow.append(.zero)
        }
      }
      identityMatrix.append(identityRow)
      diagonalMatrix.append(diagonalRow)
    }
    
    // Create column-major matrices of eigenvectors.
    let volumeElement = 1 / Double(numPoints)
    let diagonalOp = diagonalMatrix
    var phi: [[Complex<Double>]] = []
    var rhs: [[Complex<Double>]] = []
    
    for electronID in 0..<numVectors {
      var phiVector: [Complex<Double>] = []
      var rhsVector: [Complex<Double>] = []
      for cellID in 0..<numPoints {
        let expr1 = Double(electronID * cellID) * 0.1
        let expr2 = Double(cellID * (electronID + 1)) / 2
        phiVector.append(Complex.exp(Complex(0, expr1)))
        rhsVector.append(Complex(Double.cos(expr2)))
      }
      phi.append(phiVector)
      rhs.append(rhsVector)
    }
    
    func overlapDiagonal(
      _ lhs: [[Complex<Double>]],
      _ rhs: [[Complex<Double>]]
    ) -> [Complex<Double>] {
      var output: [Complex<Double>] = []
      precondition(lhs.count == 12)
      for electronID in 0..<lhs.count {
        let lhsElectron = lhs[electronID]
        let rhsElectron = rhs[electronID]
        precondition(lhsElectron.count == 100)
        var sum: Complex<Double> = .zero
        for cellID in 0..<lhsElectron.count {
          let lhsValue = lhsElectron[cellID]
          let rhsValue = rhsElectron[cellID]
          sum += lhsValue.conjugate * rhsValue * Complex(volumeElement)
        }
        output.append(sum)
      }
      return output
    }
    
    // Perform an NT matrix multiplication. Assume the left-hand side is
    // diagonal to reduce execution time in debug mode.
    func gemmDiagonal(
      _ lhs: [[Complex<Double>]],
      _ rhs: [[Complex<Double>]]
    ) -> [[Complex<Double>]] {
      var output: [[Complex<Double>]] = []
      precondition(lhs.count == 100)
      precondition(rhs.count == 12)
      for electronID in 0..<rhs.count {
        let electron = rhs[electronID]
        precondition(electron.count == 100)
        var outputVector: [Complex<Double>] = []
        for cellID in 0..<electron.count {
          let lhsElement = lhs[cellID][cellID]
          let rhsElement = electron[cellID]
          outputVector.append(lhsElement * rhsElement)
        }
        output.append(outputVector)
      }
      return output
    }
    
    func steepestDescent(
      hamiltonian: [[Complex<Double>]],
      preconditioner: [[Complex<Double>]],
      rhs: [[Complex<Double>]],
      phi: inout [[Complex<Double>]]
    ) {
      for stepID in 0..<5 {
        let HΨ = gemmDiagonal(hamiltonian, phi)
        let EΨ = rhs
        var residual = HΨ
        for electronID in 0..<12 {
          for cellID in 0..<100 {
            let lhs = HΨ[electronID][cellID]
            let rhs = EΨ[electronID][cellID]
            residual[electronID][cellID] = lhs - rhs
          }
        }
        
        let r = gemmDiagonal(preconditioner, residual)
        let Hr = gemmDiagonal(hamiltonian, r)
        let HrHr = overlapDiagonal(Hr, Hr)
        let rHr = overlapDiagonal(r, Hr)
        let rr = overlapDiagonal(r, r)
        
        // Display the norm residual.
        var output: String = ""
        output += "\(stepID)"
        precondition(phi.count == 12)
        for electronID in 0..<phi.count {
          output += "\t"
          let value = rr[electronID].real
          output += String(format: "%.3f", value)
        }
        if EigensolverTests.console {
          print(output)
        }
        
        var lambda: [Complex<Double>] = []
        for electronID in 0..<phi.count {
          let ca = HrHr[electronID]
          let cb  = 4 * rHr[electronID].real
          let cc = rr[electronID]
          
          let sqarg = Complex(cb * cb) - 4 * ca * cc
          let signB = (Complex(cb) * sqarg).real >= 0 ? Double(1) : -1
          var qq = Complex(signB) * Complex.sqrt(sqarg)
          qq = Complex(cb) + qq
          qq *= Complex(-0.5)
          lambda.append(cc / qq)
        }
        
        for electronID in 0..<12 {
          for cellID in 0..<100 {
            let lhs = phi[electronID][cellID]
            let rhs = residual[electronID][cellID]
            phi[electronID][cellID] = lhs + lambda[electronID] * rhs
          }
        }
      }
    }
    
    let numIterations = 50
    for iterationID in 0...numIterations {
      if iterationID > 0 {
        steepestDescent(
          hamiltonian: diagonalOp, preconditioner: identityMatrix,
          rhs: rhs, phi: &phi)
      }
      
      let HΨ = gemmDiagonal(diagonalOp, phi)
      let EΨ = rhs
      var residual = HΨ
      for electronID in 0..<12 {
        for cellID in 0..<100 {
          let lhs = HΨ[electronID][cellID]
          let rhs = EΨ[electronID][cellID]
          residual[electronID][cellID] = lhs - rhs
        }
      }
      let normResidual = overlapDiagonal(residual, residual)
      
      if EigensolverTests.console {
        print("  Iteration \(iterationID):")
      }
      for electronID in 0..<numVectors {
        var output = "    state \(electronID)"
        output += "  res = \(normResidual[electronID])"
        if EigensolverTests.console {
          print(output)
        }
        
        if iterationID == numIterations {
          let normres = normResidual[electronID].real
          XCTAssertLessThan(normres.magnitude, 1e-8)
        }
      }
    }
  }
  
  // Reproduce the test from the longer file:
  // https://gitlab.com/npneq/inq/-/blob/master/src/eigensolvers/steepest_descent.hpp
  //
  // Compute everything in FP32. For this test case, there is no benefit to
  // accumulating sums in FP64.
  func testDiagonalMatrixFloat() throws {
    typealias Real = Float
    
    // Create row-major hamiltonian matrices.
    let numPoints = 100
    let numVectors = 12
    var identityMatrix: [[Complex<Real>]] = []
    var diagonalMatrix: [[Complex<Real>]] = []
    
    for i in 0..<numPoints {
      var identityRow: [Complex<Real>] = []
      var diagonalRow: [Complex<Real>] = []
      for j in 0..<numPoints {
        if i == j {
          identityRow.append(Complex(1))
          diagonalRow.append(Complex(Real(i) + 1))
        } else {
          identityRow.append(.zero)
          diagonalRow.append(.zero)
        }
      }
      identityMatrix.append(identityRow)
      diagonalMatrix.append(diagonalRow)
    }
    
    // Create column-major matrices of eigenvectors.
    let volumeElement = 1 / Real(numPoints)
    let diagonalOp = diagonalMatrix
    var phi: [[Complex<Real>]] = []
    
    for electronID in 0..<numVectors {
      var phiVector: [Complex<Real>] = []
      for cellID in 0..<numPoints {
        let expr1 = Real(electronID * cellID) * 0.1
        phiVector.append(Complex.exp(Complex(0, expr1)))
      }
      phi.append(phiVector)
    }
    
    func overlapDiagonal(
      _ lhs: [[Complex<Real>]],
      _ rhs: [[Complex<Real>]]
    ) -> [Complex<Real>] {
      var output: [Complex<Real>] = []
      precondition(lhs.count == 12)
      for electronID in 0..<lhs.count {
        let lhsElectron = lhs[electronID]
        let rhsElectron = rhs[electronID]
        precondition(lhsElectron.count == 100)
        var sum: Complex<Real> = .zero
        for cellID in 0..<lhsElectron.count {
          let lhsValue = lhsElectron[cellID]
          let rhsValue = rhsElectron[cellID]
          let outValue = lhsValue.conjugate * rhsValue * Complex(volumeElement)
          sum += Complex<Real>(outValue)
        }
        output.append(Complex<Real>(sum))
      }
      return output
    }
    
    // Perform an NT matrix multiplication. Assume the left-hand side is
    // diagonal to reduce execution time in debug mode.
    func gemmDiagonal(
      _ lhs: [[Complex<Real>]],
      _ rhs: [[Complex<Real>]]
    ) -> [[Complex<Real>]] {
      var output: [[Complex<Real>]] = []
      precondition(lhs.count == 100)
      precondition(rhs.count == 12)
      for electronID in 0..<rhs.count {
        let electron = rhs[electronID]
        precondition(electron.count == 100)
        var outputVector: [Complex<Real>] = []
        for cellID in 0..<electron.count {
          let lhsElement = lhs[cellID][cellID]
          let rhsElement = electron[cellID]
          outputVector.append(lhsElement * rhsElement)
        }
        output.append(outputVector)
      }
      return output
    }
    
    func steepestDescent(
      hamiltonian: [[Complex<Real>]],
      preconditioner: [[Complex<Real>]],
      phi: inout [[Complex<Real>]]
    ) {
      var HΨ = gemmDiagonal(hamiltonian, phi)
      
      for _ in 0..<5 {
        let E = overlapDiagonal(phi, HΨ)
        let norm = overlapDiagonal(phi, phi)
        var eigenvalueNorm = E
        
        var r = HΨ
        for electronID in 0..<12 {
          eigenvalueNorm[electronID] /= norm[electronID]
          for cellID in 0..<100 {
            let E = eigenvalueNorm[electronID]
            let Ψ = phi[electronID][cellID]
            r[electronID][cellID] -= E * Ψ
          }
        }
        r = gemmDiagonal(preconditioner, r)
        
        let Hr = gemmDiagonal(hamiltonian, r)
        let rr = overlapDiagonal(r, r)
        let Ψr = overlapDiagonal(phi, r)
        let rHr = overlapDiagonal(r, Hr)
        let ΨHr = overlapDiagonal(phi, Hr)
        
        var λ: [Real] = []
        for electronID in 0..<phi.count {
          let m0 = rr[electronID]
          let m1 = Ψr[electronID]
          let m2 = rHr[electronID]
          let m3 = ΨHr[electronID]
          let m4 = E[electronID]
          let m5 = norm[electronID]
          
          let ca = (m0 * m3 - m2 * m1).real
          let cb = (m5 * m2 - m4 * m0).real
          let cc = (m4 * m1 - m3 * m5).real
          let determinant = cb + (cb * cb - 4 * ca * cc).squareRoot()
          λ.append(2 * cc / determinant)
        }
        
        for electronID in 0..<12 {
          for cellID in 0..<100 {
            let lambda = Complex(λ[electronID])
            phi[electronID][cellID] += lambda * r[electronID][cellID]
            HΨ[electronID][cellID] += lambda * Hr[electronID][cellID]
          }
        }
      }
      
      // We cannot orthogonalize here because we don't yet have a validated
      // orthogonalization operation. We also want to try an innovative approach
      // that breaks the serial dependencies in the Gram-Schmidt process. It
      // would be cross-contamination of unit tests to perform orthogonalization
      // here. Just assert that each eigenvector converges to one of the
      // expected eigenvalues.
      let norm = overlapDiagonal(phi, phi)
      for electronID in 0..<12 {
        let normalizationFactor = 1 / norm[electronID].real.squareRoot()
        for cellID in 0..<100 {
          phi[electronID][cellID] *= Complex(normalizationFactor)
        }
      }
    }
    
    let numIterations = 50
    for iterationID in 0...numIterations {
      if iterationID > 0 {
        steepestDescent(
          hamiltonian: diagonalOp, preconditioner: identityMatrix, phi: &phi)
      }
      
      let HΨ = gemmDiagonal(diagonalOp, phi)
      let E = overlapDiagonal(phi, HΨ)
      var residual = HΨ
      for electronID in 0..<12 {
        for cellID in 0..<100 {
          let Ψ = phi[electronID][cellID]
          residual[electronID][cellID] -= E[electronID] * Ψ
        }
      }
      let normResidual = overlapDiagonal(residual, residual)
      
      if EigensolverTests.console {
        print("  Iteration \(iterationID):")
      }
      for electronID in 0..<numVectors {
        var output = "    state \(electronID)"
        output += "  evalue = \(E[electronID])"
        output += "  res = \(normResidual[electronID])"
        if EigensolverTests.console {
          print(output)
        }
        
        if iterationID == numIterations {
          let normres = normResidual[electronID].real
          XCTAssertLessThan(normres.magnitude, 1e-3)
        }
      }
    }
  }
}
