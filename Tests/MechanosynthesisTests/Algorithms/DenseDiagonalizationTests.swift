import XCTest
import Accelerate // Gate out once there are Swift kernels for CPUs without AMX.
import Numerics
import QuartzCore

// TODO: Create a custom direct diagonalization kernel. Start with the
// simplest one known, then move on to more advanced solvers.
final class DenseDiagonalizationTests: XCTestCase {
  // Standalone function for tridiagonalizing a matrix.
  // - Performs several loop iterations with O(n^2) complexity each. The overall
  //   compute cost is O(n^3).
  static func tridiagonalize(
    matrix originalMatrix: [Float],
    n: Int
  ) -> [Float] {
    var currentMatrixA = originalMatrix
    
    // This requires that n > 1.
    for k in 0..<n - 2 {
      // Determine 'α' and 'r'.
      let address_k1k = (k + 1) * n + k
      let ak1k = currentMatrixA[address_k1k]
      var alpha: Float = .zero
      
      for rowID in (k + 1)..<n {
        let columnID = k
        let address = rowID * n + columnID
        let value = currentMatrixA[address]
        alpha += value * value
      }
      alpha.formSquareRoot()
      alpha *= (ak1k >= 0) ? -1 : 1
      
      var r = alpha * alpha - ak1k * alpha
      r = (r / 2).squareRoot()
      
      // Construct 'v'.
      var v = [Float](repeating: 0, count: n)
      v[k + 1] = (ak1k - alpha) / (2 * r)
      for vectorLane in (k + 2)..<n {
        let matrixAddress = vectorLane * n + k
        let matrixValue = currentMatrixA[matrixAddress]
        v[vectorLane] = matrixValue / (2 * r)
      }
      
      // Operation 1: gemv(A, v)
      var operationResult1 = [Float](repeating: 0, count: n)
      for rowID in 0..<n {
        var dotProduct: Float = .zero
        for columnID in 0..<n {
          let matrixAddress = rowID * n + columnID
          let vectorAddress = columnID
          let matrixValue = currentMatrixA[matrixAddress]
          let vectorValue = v[vectorAddress]
          dotProduct += matrixValue * vectorValue
        }
        operationResult1[rowID] = dotProduct
      }
      
      // Operation 2: scatter(..., v)
      for rowID in 0..<n {
        for columnID in 0..<n {
          let rowValue = operationResult1[rowID]
          let columnValue = v[columnID]
          let matrixAddress = rowID * n + columnID
          currentMatrixA[matrixAddress] -= 2 * rowValue * columnValue
        }
      }
      
      // Operation 3: gemv(vT, A)
      var operationResult3 = [Float](repeating: 0, count: n)
      for columnID in 0..<n {
        var dotProduct: Float = .zero
        for rowID in 0..<n {
          let vectorAddress = rowID
          let matrixAddress = rowID * n + columnID
          let vectorValue = v[vectorAddress]
          let matrixValue = currentMatrixA[matrixAddress]
          dotProduct += vectorValue * matrixValue
        }
        operationResult3[columnID] = dotProduct
      }
      
      // Operation 4: scatter(v, ...)
      for rowID in 0..<n {
        for columnID in 0..<n {
          let rowValue = v[rowID]
          let columnValue = operationResult3[columnID]
          let matrixAddress = rowID * n + columnID
          currentMatrixA[matrixAddress] -= 2 * rowValue * columnValue
        }
      }
    }
    return currentMatrixA
  }
}
