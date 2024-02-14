import XCTest
import DFT
import Numerics

// Test the steepest descent algorithm from INQ.
final class EigensolverTests: XCTestCase {
  // Reproduce the test from the shorter file:
  // https://gitlab.com/npneq/inq/-/blob/master/src/solvers/steepest_descent.hpp
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
        rhsVector.append(Complex(Double.cos(expr1)))
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
    
    
  }
  
  // Reproduce the test from the longer file:
  // https://gitlab.com/npneq/inq/-/blob/master/src/eigensolvers/steepest_descent.hpp
  func testDiagonalMatrixFloat() throws {
    
  }
}
