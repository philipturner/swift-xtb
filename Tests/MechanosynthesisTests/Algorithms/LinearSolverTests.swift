import XCTest
import Mechanosynthesis
import Numerics

final class LinearSolverTests: XCTestCase {
  static let gridSize: Int = 4
  static let h: Float = 0.5
  static var cellCount: Int { gridSize * gridSize * gridSize }
  
  // Create the 'b' vector, which equals -4πρ.
  static func createScaledChargeDensity() -> [Float] {
    var output = [Float](repeating: .zero, count: cellCount)
    for permutationZ in -1...0 {
      for permutationY in -1...0 {
        for permutationX in -1...0 {
          var indices = SIMD3(repeating: gridSize / 2)
          indices[0] += permutationX
          indices[1] += permutationY
          indices[2] += permutationZ
          
          // Place 1/8 of the charge density in each of the 8 cells.
          let chargeEnclosed: Float = 1.0 / 8
          let chargeDensity = chargeEnclosed / (h * h * h)
          
          // Multiply -4π with ρ, resulting in -4πρ.
          let rhsValue = (-4 * Float.pi) * chargeDensity
          
          // Write the right-hand side to memory.
          var cellID = indices.z * (gridSize * gridSize)
          cellID += indices.y * gridSize + indices.x
          output[cellID] = rhsValue
        }
      }
    }
    return output
  }
  
  static func createAddress(indices: SIMD3<Int>) -> Int {
    indices.z * (gridSize * gridSize) + indices.y * gridSize + indices.x
  }
  
  // Apply the 'A' matrix, which equals ∇^2.
  //
  // The Laplacian has second-order accuracy. It fills ghost cells with the
  // multipole expansion of the charge enclosed.
  static func applyLaplacian(_ x: [Float]) -> [Float] {
    guard x.count == cellCount else {
      fatalError("Dimensions of 'x' did not match problem size.")
    }
    
    // Iterate over the cells.
    var output = [Float](repeating: 0, count: cellCount)
    for indexZ in 0..<gridSize {
      for indexY in 0..<gridSize {
        for indexX in 0..<gridSize {
          var dotProduct: Float = .zero
          
          // Apply the FMA on the diagonal.
          let cellIndices = SIMD3(indexX, indexY, indexZ)
          let cellAddress = createAddress(indices: cellIndices)
          let cellValue = x[cellAddress]
          dotProduct += -6 / (h * h) * cellValue
          
          // Iterate over the faces.
          for faceID in 0..<6 {
            let coordinateID = faceID / 2
            let coordinateShift = (faceID % 2 == 0) ? -1 : 1
            
            // Locate the neighboring cell.
            var neighborIndices = SIMD3(indexX, indexY, indexZ)
            neighborIndices[coordinateID] += coordinateShift
            
            if all(neighborIndices .>= 0) && all(neighborIndices .< gridSize) {
              let neighborAddress = createAddress(indices: neighborIndices)
              let neighborValue = x[neighborAddress]
              dotProduct += 1 / (h * h) * neighborValue
            } else {
              var neighborPosition = SIMD3<Float>(neighborIndices)
              neighborPosition = h * (neighborPosition + 0.5)
              var nucleusPosition = SIMD3(repeating: Float(gridSize))
              nucleusPosition = h * (nucleusPosition * 0.5)
              
              // Generate a ghost value from the point charge approximation.
              let r = neighborPosition - nucleusPosition
              let distance = (r * r).sum().squareRoot()
              let neighborValue = 1 / distance
              dotProduct += 1 / (h * h) * neighborValue
            }
          }
          
          // Store the dot product.
          output[cellAddress] = dotProduct
        }
      }
    }
    
    return output
  }
  
  // Create the analytical value for the solution.
  static func createReferenceSolution() -> [Float] {
    var output = [Float](repeating: .zero, count: Self.cellCount)
    for indexZ in 0..<Self.gridSize {
      for indexY in 0..<Self.gridSize {
        for indexX in 0..<Self.gridSize {
          let cellIndices = SIMD3(indexX, indexY, indexZ)
          let cellAddress = Self.createAddress(indices: cellIndices)
          
          var cellPosition = SIMD3<Float>(cellIndices)
          cellPosition = Self.h * (cellPosition + 0.5)
          var nucleusPosition = SIMD3(repeating: Float(Self.gridSize))
          nucleusPosition = Self.h * (nucleusPosition * 0.5)
          
          // Generate a ghost value from the point charge approximation.
          let r = cellPosition - nucleusPosition
          let distance = (r * r).sum().squareRoot()
          let cellValue = 1 / distance
          
          // Store the dot product.
          output[cellAddress] = cellValue
        }
      }
    }
    return output
  }
  
  // Shift a vector by a constant times another vector.
  //
  // Returns: original + scale * correction
  static func shift(
    _ original: [Float],
    scale: Float,
    correction: [Float]
  ) -> [Float] {
    var output = [Float](repeating: .zero, count: Self.cellCount)
    for cellID in 0..<Self.cellCount {
      var cellValue = original[cellID]
      cellValue += scale * correction[cellID]
      output[cellID] = cellValue
    }
    return output
  }
  
  // Take the dot product of two vectors.
  static func dot(
    _ lhs: [Float],
    _ rhs: [Float]
  ) -> Float {
    var accumulator: Double = .zero
    for cellID in 0..<Self.cellCount {
      let lhsValue = lhs[cellID]
      let rhsValue = rhs[cellID]
      accumulator += Double(lhsValue * rhsValue)
    }
    return Float(accumulator)
  }
  
  // Jacobi method:
  //
  // Ax = b
  // (D + L + U)x = b
  // Dx = b - (L + U)x
  // Dx = b - (A - D)x
  // Dx = b - Ax + Dx
  // x = x + D^{-1} (b - Ax)
  func testJacobiMethod() throws {
    let b = Self.createScaledChargeDensity()
    var x = [Float](repeating: .zero, count: Self.cellCount)
    let expectedX = Self.createReferenceSolution()
    
    for _ in 0..<20 {
      let Ax = Self.applyLaplacian(x)
      let r = Self.shift(b, scale: -1, correction: Ax)
      
      do {
        let r2 = Self.dot(r, r)
        let resNorm = r2.squareRoot()
        print("||r|| = \(resNorm)")
      }
      
      let D = -6 / (Self.h * Self.h)
      x = Self.shift(x, scale: 1 / D, correction: r)
    }
  }
}
