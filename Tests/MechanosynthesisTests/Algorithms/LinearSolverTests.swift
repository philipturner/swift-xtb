import XCTest
import Mechanosynthesis
import Numerics

final class LinearSolverTests: XCTestCase {
  static let gridSize: Int = 8
  static let h: Float = 0.25
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
  
  // Apply the 'A' matrix (∇^2), while omitting ghost cells.
  //
  // The Laplacian has second-order accuracy.
  static func applyLaplacianLinearPart(_ x: [Float]) -> [Float] {
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
              // dotProduct += 1 / (h * h) * neighborValue
            }
          }
          
          // Store the dot product.
          output[cellAddress] = dotProduct
        }
      }
    }
    
    return output
  }
  
  // The Laplacian, omitting contributions from the input 'x'.
  //
  // Fills ghost cells with the multipole expansion of the charge enclosed.
  static func applyLaplacianBoundary() -> [Float] {
    // Iterate over the cells.
    var output = [Float](repeating: 0, count: cellCount)
    for indexZ in 0..<gridSize {
      for indexY in 0..<gridSize {
        for indexX in 0..<gridSize {
          var dotProduct: Float = .zero
          
          // Apply the FMA on the diagonal.
          let cellIndices = SIMD3(indexX, indexY, indexZ)
          let cellAddress = createAddress(indices: cellIndices)
          
          // Iterate over the faces.
          for faceID in 0..<6 {
            let coordinateID = faceID / 2
            let coordinateShift = (faceID % 2 == 0) ? -1 : 1
            
            // Locate the neighboring cell.
            var neighborIndices = SIMD3(indexX, indexY, indexZ)
            neighborIndices[coordinateID] += coordinateShift
            
            if all(neighborIndices .>= 0) && all(neighborIndices .< gridSize) {
              let neighborAddress = createAddress(indices: neighborIndices)
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
  
  // MARK: - Utilities
  
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
  
  // MARK: - Tests
  
  // Jacobi method:
  //
  // Ax = b
  // (D + L + U)x = b
  // Dx = b - (L + U)x
  // Dx = b - (A - D)x
  // Dx = b - Ax + Dx
  // x = x + D^{-1} (b - Ax)
  func testJacobiMethod() throws {
    var b = Self.createScaledChargeDensity()
    let L2x = Self.applyLaplacianBoundary()
    b = Self.shift(b, scale: -1, correction: L2x)
    
    var x = [Float](repeating: .zero, count: Self.cellCount)
    print()
    print("Jacobi")
    for _ in 0..<30 {
      let L1x = Self.applyLaplacianLinearPart(x)
      let r = Self.shift(b, scale: -1, correction: L1x)
      let r2 = Self.dot(r, r)
      let resNorm = r2.squareRoot()
      print("||r|| = \(resNorm)")
      
      let D = -6 / (Self.h * Self.h)
      x = Self.shift(x, scale: 1 / D, correction: r)
    }
  }
  
  // Conjugate gradient method:
  //
  // r = b - Ax
  // p = r - Σ_i < p_i | A | r > / < p_i | A | p_i >
  // a = < p | r > / < p | A | p >
  // x = x + a p
  //
  // Efficient version:
  //
  // r = b - Ax
  // p = r
  // repeat
  //   a = < r | r > / < p | A | p >
  //   x_new = x + a p
  //   r_new = r - a A p
  //
  //   b = < r_new | r_new > / < r | r >
  //   p_new = r_new + b p
  func testConjugateGradientMethod() throws {
    var b = Self.createScaledChargeDensity()
    let L2x = Self.applyLaplacianBoundary()
    b = Self.shift(b, scale: -1, correction: L2x)
    
    func logNormRes() {
      let L1x = Self.applyLaplacianLinearPart(x)
      let r = Self.shift(b, scale: -1, correction: L1x)
      let r2 = Self.dot(r, r)
      let resNorm = r2.squareRoot()
      print("||r|| = \(resNorm)")
    }
    
    var x = [Float](repeating: .zero, count: Self.cellCount)
    let L1x = Self.applyLaplacianLinearPart(x)
    var r = Self.shift(b, scale: -1, correction: L1x)
    var p = r
    var rr = Self.dot(r, r)
    print()
    print("Conjugate Gradient")
    logNormRes()
    
    for _ in 0..<30 {
      let Ap = Self.applyLaplacianLinearPart(p)
      
      let a = rr / Self.dot(p, Ap)
      let xNew = Self.shift(x, scale: a, correction: p)
      let rNew = Self.shift(r, scale: -a, correction: Ap)
      let rrNew = Self.dot(rNew, rNew)
      
      let b = rrNew / rr
      let pNew = Self.shift(rNew, scale: b, correction: p)
      
      x = xNew
      r = rNew
      p = pNew
      rr = rrNew
      
      logNormRes()
    }
  }
  
  // Gauss-Seidel method:
  //
  // x_i = (1 / a_ii) (b_i - Σ_(j ≠ i) a_ij x_j)
  //
  // Red-black scheme:
  //
  // iterate over all the red cells in parallel
  // iterate over all the black cells in parallel
  // only works with 2nd order FD
  func testGaussSeidelMethod() {
    var b = Self.createScaledChargeDensity()
    let L2x = Self.applyLaplacianBoundary()
    b = Self.shift(b, scale: -1, correction: L2x)
    
    print()
    print("Gauss-Seidel")
    var x = [Float](repeating: .zero, count: Self.cellCount)
    
    // Updates all of the selected cells in-place.
    func executeSweep(red: Bool, black: Bool) {
      for indexZ in 0..<Self.gridSize {
        for indexY in 0..<Self.gridSize {
          for indexX in 0..<Self.gridSize {
            var dotProduct: Float = .zero
            
            // Mask out either the red or black cells.
            let parity = indexX ^ indexY ^ indexZ
            switch parity & 1 {
            case 0:
              guard red else {
                continue
              }
            case 1:
              guard black else {
                continue
              }
            default:
              fatalError("This should never happen.")
            }
            
            // Iterate over the faces.
            for faceID in 0..<6 {
              let coordinateID = faceID / 2
              let coordinateShift = (faceID % 2 == 0) ? -1 : 1
              
              // Locate the neighboring cell.
              var neighborIndices = SIMD3(indexX, indexY, indexZ)
              neighborIndices[coordinateID] += coordinateShift
              
              if all(neighborIndices .>= 0),
                 all(neighborIndices .< Self.gridSize) {
                let neighborAddress = Self
                  .createAddress(indices: neighborIndices)
                let neighborValue = x[neighborAddress]
                dotProduct += 1 / (Self.h * Self.h) * neighborValue
              }
            }
            
            let cellIndices = SIMD3(indexX, indexY, indexZ)
            let cellAddress = Self.createAddress(indices: cellIndices)
            
            // Overwrite the current value.
            let rhsValue = b[cellAddress]
            let diagonalValue: Float = -6 / (Self.h * Self.h)
            let newValue = (rhsValue - dotProduct) / diagonalValue
            x[cellAddress] = newValue
          }
        }
      }
    }
    
    for iterationID in 0..<30 {
      let L1x = Self.applyLaplacianLinearPart(x)
      let r = Self.shift(b, scale: -1, correction: L1x)
      let r2 = Self.dot(r, r)
      let resNorm = r2.squareRoot()
      print("||r|| = \(resNorm)")
      
      executeSweep(red: true, black: false)
      executeSweep(red: false, black: true)
    }
  }
  
  func testMultigridMethod() {
    print()
    print("Multigrid")
    var b = Self.createScaledChargeDensity()
    var x = [Float](repeating: .zero, count: Self.cellCount)
    
    func logNormRes() {
      let L1x = Self.applyLaplacianLinearPart(x)
      let L2x = Self.applyLaplacianBoundary()
      let Ax = Self.shift(L1x, scale: 1, correction: L2x)
      var r = Self.shift(b, scale: -1, correction: Ax)
      do {
        let r2 = Self.dot(r, r)
        let resNorm = r2.squareRoot()
        print("||r|| = \(resNorm)")
      }
    }
    
    for iterationID in 0..<30 {
      logNormRes()
      
      
      // Gauss-Seidel red-black
      func GSRB_LEVEL(
        e: inout [Float], r: [Float], red: Bool
      ) {
        for indexZ in 0..<Self.gridSize {
          for indexY in 0..<Self.gridSize {
            for indexX in 0..<Self.gridSize {
              // Mask out either the red or black cells.
              let parity = indexX ^ indexY ^ indexZ
              switch parity & 1 {
              case 0:
                guard red else {
                  continue
                }
              case 1:
                guard !red else {
                  continue
                }
              default:
                fatalError("This should never happen.")
              }
              
              let cellIndices = SIMD3(indexX, indexY, indexZ)
              let cellAddress = Self.createAddress(indices: cellIndices)
              
              // Iterate over the faces.
              var faceAccumulator: Float = .zero
              for faceID in 0..<6 {
                let coordinateID = faceID / 2
                let coordinateShift = (faceID % 2 == 0) ? -1 : 1
                
                // Locate the neighboring cell.
                var neighborIndices = SIMD3(indexX, indexY, indexZ)
                neighborIndices[coordinateID] += coordinateShift
                guard all(neighborIndices .>= 0),
                      all(neighborIndices .< Self.gridSize) else {
                  // Add 'zero' to the accumulator.
                  continue
                }
                
                // Add the neighbor's value to the accumulator.
                let neighborAddress = Self
                  .createAddress(indices: neighborIndices)
                let neighborValue = e[neighborAddress]
                faceAccumulator += 1 / (Self.h * Self.h) * neighborValue
              }
              
              // Fetch the values to evaluate GSRB_LEVEL(e, R, h).
              let rValue = r[cellAddress]
              var eValue = e[cellAddress]
              
              // Update the error in-place.
              var λ = Self.h * Self.h / 6
              e[cellAddress] = λ * (faceAccumulator - rValue)
            }
          }
        }
      }
      
      // Emulate a multigrid V-cycle, except the "coarse" levels are actually
      // as fine as the "fine" levels.
      var rFine: [Float]
      do {
        let L1x = Self.applyLaplacianLinearPart(x)
        let L2x = Self.applyLaplacianBoundary()
        let Ax = Self.shift(L1x, scale: 1, correction: L2x)
        rFine = Self.shift(b, scale: -1, correction: Ax)
      }
      
      // Smoothing iterations on the fine level.
      var eFine = [Float](repeating: .zero, count: Self.cellCount)
      GSRB_LEVEL(e: &eFine, r: rFine, red: true)
      GSRB_LEVEL(e: &eFine, r: rFine, red: false)
      
      // Emulate the averaging from fine -> coarse.
      var rCoarse = Self.shift(
        rFine, scale: -1, correction: Self.applyLaplacianLinearPart(eFine))
      
      // Smoothing iterations on the coarse level.
      var eCoarse = [Float](repeating: .zero, count: Self.cellCount)
//      GSRB_LEVEL(e: &eCoarse, r: rCoarse, red: true)
//      GSRB_LEVEL(e: &eCoarse, r: rCoarse, red: false)
//      GSRB_LEVEL(e: &eCoarse, r: rCoarse, red: true)
//      GSRB_LEVEL(e: &eCoarse, r: rCoarse, red: false)
      
      // Emulate the correction from coarse -> fine.
      eFine = Self.shift(eFine, scale: 1, correction: eCoarse)
      rFine = Self.shift(rFine, scale: -1, correction: Self.applyLaplacianLinearPart(eFine))
      
      // Smoothing iterations on the fine level.
      var δeFine = [Float](repeating: .zero, count: Self.cellCount)
      GSRB_LEVEL(e: &δeFine, r: rFine, red: true)
      GSRB_LEVEL(e: &δeFine, r: rFine, red: false)
      
      // Add the error vector to the solution.
      eFine = Self.shift(eFine, scale: 1, correction: δeFine)
      x = Self.shift(x, scale: 1, correction: eFine)
    }
  }
}
