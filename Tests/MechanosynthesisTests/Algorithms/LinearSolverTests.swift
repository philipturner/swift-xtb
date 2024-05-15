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
    for _ in 0..<30 {
      let L1x = Self.applyLaplacianLinearPart(x)
      let r = Self.shift(b, scale: -1, correction: L1x)
      
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
    // TODO: Compare CG to MG again, once the preconditioner is working. How
    // does it perform when the grid expands from 8x8x8 to 32x32x32?
    
    var b = Self.createScaledChargeDensity()
    let L2x = Self.applyLaplacianBoundary()
    b = Self.shift(b, scale: -1, correction: L2x)
    
    var x = [Float](repeating: .zero, count: Self.cellCount)
    let L1x = Self.applyLaplacianLinearPart(x)
    var r = Self.shift(b, scale: -1, correction: L1x)
    var p = r
    var rr = Self.dot(r, r)
    
    print()
    print("Conjugate Gradient")
    for _ in 0..<30 {
      do {
        let L1x = Self.applyLaplacianLinearPart(x)
        let r = Self.shift(b, scale: -1, correction: L1x)
        let r2 = Self.dot(r, r)
        let normres = r2.squareRoot()
        print("||r|| = \(normres)")
      }
      
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
    }
  }
  
  // Preconditioned conjugate gradient method:
  //
  // r = b - Ax
  // p = K r
  // repeat
  //   a = < r | K | r > / < p | A | p >
  //   x_new = x + a p
  //   r_new = r - a A p
  //
  //   b = < r_new | K | r_new > / < r | K | r >
  //   p_new = K r_new + b p
  func testPreconditionedConjugateGradient() {
    // First, testing steepest descent as described in the real-space DFT
    // textbook. This seems slightly different than Jacobi iteration, which
    // explicitly fetches the diagonal entries and inverts them.
    var b = Self.createScaledChargeDensity()
    let L2x = Self.applyLaplacianBoundary()
    b = Self.shift(b, scale: -1, correction: L2x)
    
    var x = [Float](repeating: .zero, count: Self.cellCount)
    let L1x = Self.applyLaplacianLinearPart(x)
    var r = Self.shift(b, scale: -1, correction: L1x)
    var p = applyLaplacianPreconditioner(r)
    
    print()
    print("Preconditioned Conjugate Gradient")
    for _ in 0..<15 {
      do {
        let L1x = Self.applyLaplacianLinearPart(x)
        let r = Self.shift(b, scale: -1, correction: L1x)
        let r2 = Self.dot(r, r)
        let normres = r2.squareRoot()
        print("||r|| = \(normres)")
      }
      
      let a = Self.dot(r, applyLaplacianPreconditioner(r)) / Self.dot(p, Self.applyLaplacianLinearPart(p))
      let xNew = Self.shift(x, scale: a, correction: p)
      let rNew = Self.shift(r, scale: -a, correction: Self.applyLaplacianLinearPart(p))
      
      let b = Self.dot(rNew, applyLaplacianPreconditioner(rNew)) / Self.dot(r, applyLaplacianPreconditioner(r))
      let pNew = Self.shift(applyLaplacianPreconditioner(rNew), scale: b, correction: p)
      
      x = xNew
      r = rNew
      p = pNew
    }
    
    func applyLaplacianPreconditioner(_ x: [Float]) -> [Float] {
      let h = Self.h
      let gridSize = Self.gridSize
      let cellCount = Self.cellCount
      
      @_transparent
      func createAddress(indices: SIMD3<Int16>) -> Int {
        Int(indices.z) * (gridSize * gridSize) +
        Int(indices.y) * gridSize +
        Int(indices.x)
      }
      
      // Pre-compile a list of neighbor offsets.
      var neighborData: [SIMD4<Int16>] = []
      for offsetZ in -2...2 {
        for offsetY in -2...2 {
          for offsetX in -2...2 {
            var indices = SIMD3(Int16(offsetX), Int16(offsetY), Int16(offsetZ))
            let integerDistanceSquared = (indices &* indices).wrappedSum()
            
            // This tolerance creates a 33-point convolution kernel.
            guard integerDistanceSquared <= 4 else {
              continue
            }
            
            // Execute the formula for matrix elements.
            var K: Float = .zero
            K += 0.6 * Float.exp(-2.25 * Float(integerDistanceSquared))
            K += 0.4 * Float.exp(-0.72 * Float(integerDistanceSquared))
            let quantized = Int16(K * 32767)
            
            // Pack the data into a compact 64-bit word.
            let vector = SIMD4(indices, quantized)
            neighborData.append(vector)
          }
        }
      }
      
      // Iterate over the cells.
      var output = [Float](repeating: 0, count: cellCount)
      for indexZ in 0..<gridSize {
        for indexY in 0..<gridSize {
          for indexX in 0..<gridSize {
            // Iterate over the convolution points.
            var accumulator: Float = .zero
            
            // The test took 0.015 seconds before.
            // 0.013 seconds
            let cellIndices64 = SIMD3(indexX, indexY, indexZ)
            let cellIndices = SIMD3<Int16>(truncatingIfNeeded: cellIndices64)
            for vector in neighborData {
              let offset = unsafeBitCast(vector, to: SIMD3<Int16>.self)
              let neighborIndices = cellIndices &+ offset
              guard all(neighborIndices .>= 0),
                    all(neighborIndices .< Int16(gridSize)) else {
                continue
              }
              
              // Read the neighbor data point from memory.
              let neighborAddress = createAddress(indices: neighborIndices)
              let neighborValue = x[neighborAddress]
              let K = Float(vector[3]) / 32767
              accumulator += neighborValue * K
            }
            
            // Write the convolution result to memory.
            let cellAddress = createAddress(indices: cellIndices)
            output[cellAddress] = accumulator
          }
        }
      }
      
      return output
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
      executeSweep(red: true, black: false)
      executeSweep(red: false, black: true)
    }
  }
  
  func testMultigridMethod() {
    var b = Self.createScaledChargeDensity()
    var x = [Float](repeating: .zero, count: Self.cellCount)
    
    // One V-cycle should be treated as two SD or CG iterations.
    for iterationID in 0..<15 {
      // Initialize the residual.
      var rFine: [Float]
      do {
        let L1x = Self.applyLaplacianLinearPart(x)
        let L2x = Self.applyLaplacianBoundary()
        let Ax = Self.shift(L1x, scale: 1, correction: L2x)
        rFine = Self.shift(b, scale: -1, correction: Ax)
      }
      
      // Smoothing iterations on the fine level.
      var eFine = [Float](repeating: .zero, count: Self.cellCount)
      GSRB_LEVEL(e: &eFine, r: rFine, coarseness: 1, red: true)
      GSRB_LEVEL(e: &eFine, r: rFine, coarseness: 1, red: false)
      
      // Encapsulate the coarse level.
      if true {
        // Restrict from fine to coarse.
        var rFineCorrected = rFine
        var rCoarse: [Float] = []
        correctResidual(e: eFine, r: &rFineCorrected, coarseness: 1)
        shiftResolution(
          fineGrid: &rFineCorrected, coarseGrid: &rCoarse,
          fineLevelCoarseness: 1, shiftingUp: true)
        
        // Smoothing iterations on the coarse level.
        var eCoarse = [Float](repeating: .zero, count: Self.cellCount / 8)
        GSRB_LEVEL(e: &eCoarse, r: rCoarse, coarseness: 2, red: true)
        GSRB_LEVEL(e: &eCoarse, r: rCoarse, coarseness: 2, red: false)
        
        // Encapsulate the tiny level.
        if true {
          // Restrict from coarse to tiny.
          var rCoarseCorrected = rCoarse
          var rTiny: [Float] = []
          correctResidual(e: eCoarse, r: &rCoarseCorrected, coarseness: 2)
          shiftResolution(
            fineGrid: &rCoarseCorrected, coarseGrid: &rTiny,
            fineLevelCoarseness: 2, shiftingUp: true)
          
          // Smoothing iterations on the tiny level.
          //
          // NOTE: When the tiny grid is 2x2x2, smoothing iterations on this
          // level actually harm convergence. It only helps when the tiny grid
          // is 4x4x4 or larger. However, the test covers this level as a proof
          // of concept.
          var eTiny = [Float](repeating: .zero, count: Self.cellCount / 64)
          GSRB_LEVEL(e: &eTiny, r: rTiny, coarseness: 4, red: true)
          GSRB_LEVEL(e: &eTiny, r: rTiny, coarseness: 4, red: false)
          
          // Prolong from tiny to coarse.
          shiftResolution(
            fineGrid: &eCoarse, coarseGrid: &eTiny,
            fineLevelCoarseness: 2, shiftingUp: false)
          correctResidual(e: eCoarse, r: &rCoarse, coarseness: 2)
          
          // Smoothing iterations on the coarse level.
          var δeCoarse = [Float](repeating: .zero, count: Self.cellCount / 8)
          GSRB_LEVEL(e: &δeCoarse, r: rCoarse, coarseness: 2, red: true)
          GSRB_LEVEL(e: &δeCoarse, r: rCoarse, coarseness: 2, red: false)
          for cellID in 0..<Self.cellCount / 8 {
            eCoarse[cellID] += δeCoarse[cellID]
          }
        }
        
        // Prolong from coarse to fine.
        shiftResolution(
          fineGrid: &eFine, coarseGrid: &eCoarse,
          fineLevelCoarseness: 1, shiftingUp: false)
        correctResidual(e: eFine, r: &rFine, coarseness: 1)
        
        // Smoothing iterations on the fine level.
        var δeFine = [Float](repeating: .zero, count: Self.cellCount)
        GSRB_LEVEL(e: &δeFine, r: rFine, coarseness: 1, red: true)
        GSRB_LEVEL(e: &δeFine, r: rFine, coarseness: 1, red: false)
        for cellID in 0..<Self.cellCount {
          eFine[cellID] += δeFine[cellID]
        }
      }
      
      // Update the solution.
      x = Self.shift(x, scale: 1, correction: eFine)
    }
    
    // Gauss-Seidel red-black
    func GSRB_LEVEL(
      e: inout [Float], r: [Float], coarseness: Int, red: Bool
    ) {
      let h = Self.h * Float(coarseness)
      let gridSize = Self.gridSize / coarseness
      func createAddress(indices: SIMD3<Int>) -> Int {
        indices.z * (gridSize * gridSize) + indices.y * gridSize + indices.x
      }
      
      for indexZ in 0..<gridSize {
        for indexY in 0..<gridSize {
          for indexX in 0..<gridSize {
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
            let cellAddress = createAddress(indices: cellIndices)
            
            // Iterate over the faces.
            var faceAccumulator: Float = .zero
            for faceID in 0..<6 {
              let coordinateID = faceID / 2
              let coordinateShift = (faceID % 2 == 0) ? -1 : 1
              
              // Locate the neighboring cell.
              var neighborIndices = SIMD3(indexX, indexY, indexZ)
              neighborIndices[coordinateID] += coordinateShift
              guard all(neighborIndices .>= 0),
                    all(neighborIndices .< gridSize) else {
                // Add 'zero' to the accumulator.
                continue
              }
              
              // Add the neighbor's value to the accumulator.
              let neighborAddress = createAddress(indices: neighborIndices)
              let neighborValue = e[neighborAddress]
              faceAccumulator += 1 / (h * h) * neighborValue
            }
            
            // Fetch the values to evaluate GSRB_LEVEL(e, R, h).
            let rValue = r[cellAddress]
            var eValue = e[cellAddress]
            
            // Update the error in-place.
            var λ = h * h / 6
            e[cellAddress] = λ * (faceAccumulator - rValue)
          }
        }
      }
    }
    
    func correctResidual(e: [Float], r: inout [Float], coarseness: Int) {
      let h = Self.h * Float(coarseness)
      let gridSize = Self.gridSize / coarseness
      func createAddress(indices: SIMD3<Int>) -> Int {
        indices.z * (gridSize * gridSize) + indices.y * gridSize + indices.x
      }
      
      // Iterate over the cells.
      for indexZ in 0..<gridSize {
        for indexY in 0..<gridSize {
          for indexX in 0..<gridSize {
            var dotProduct: Float = .zero
            
            // Apply the FMA on the diagonal.
            let cellIndices = SIMD3(indexX, indexY, indexZ)
            let cellAddress = createAddress(indices: cellIndices)
            let cellValue = e[cellAddress]
            dotProduct += -6 / (h * h) * cellValue
            
            // Iterate over the faces.
            for faceID in 0..<6 {
              let coordinateID = faceID / 2
              let coordinateShift = (faceID % 2 == 0) ? -1 : 1
              
              // Locate the neighboring cell.
              var neighborIndices = SIMD3(indexX, indexY, indexZ)
              neighborIndices[coordinateID] += coordinateShift
              guard all(neighborIndices .>= 0),
                    all(neighborIndices .< gridSize) else {
                // Add 'zero' to the dot product.
                continue
              }
              
              let neighborAddress = createAddress(indices: neighborIndices)
              let neighborValue = e[neighborAddress]
              dotProduct += 1 / (h * h) * neighborValue
            }
            
            // Update the residual.
            let L2e = dotProduct
            r[cellAddress] -= L2e
          }
        }
      }
    }
    
    // Performs a power-2 shift to a coarser level.
    func shiftResolution(
      fineGrid: inout [Float], coarseGrid: inout [Float],
      fineLevelCoarseness: Int, shiftingUp: Bool
    ) {
      let fineGridSize = Self.gridSize / fineLevelCoarseness
      let coarseGridSize = fineGridSize / 2
      let coarseCellCount = coarseGridSize * coarseGridSize * coarseGridSize
      
      func createCoarseAddress(indices: SIMD3<Int>) -> Int {
        indices.z * (coarseGridSize * coarseGridSize) +
        indices.y * coarseGridSize + indices.x
      }
      func createFineAddress(indices: SIMD3<Int>) -> Int {
        indices.z * (fineGridSize * fineGridSize) +
        indices.y * fineGridSize + indices.x
      }
      
      // Create the coarse grid.
      if shiftingUp {
        coarseGrid = [Float](repeating: .zero, count: coarseCellCount)
      }
      
      // Iterate over the coarse grid.
      for indexZ in 0..<coarseGridSize {
        for indexY in 0..<coarseGridSize {
          for indexX in 0..<coarseGridSize {
            // Read from the coarse grid.
            let coarseIndices = SIMD3<Int>(indexX, indexY, indexZ)
            let coarseAddress = createCoarseAddress(indices: coarseIndices)
            let coarseValue = coarseGrid[coarseAddress]
            
            // Iterate over the footprint on the finer grid.
            var accumulator: Float = .zero
            for permutationZ in 0..<2 {
              for permutationY in 0..<2 {
                for permutationX in 0..<2 {
                  var fineIndices = 2 &* coarseIndices
                  fineIndices[0] += permutationX
                  fineIndices[1] += permutationY
                  fineIndices[2] += permutationZ
                  let fineAddress = createFineAddress(indices: fineIndices)
                  
                  if shiftingUp {
                    // Read from the fine grid.
                    let fineValue = fineGrid[fineAddress]
                    accumulator += (1.0 / 8) * fineValue
                  } else {
                    // Update the fine grid.
                    fineGrid[fineAddress] += coarseValue
                  }
                }
              }
            }
            
            // Update the coarse grid.
            if shiftingUp {
              coarseGrid[coarseAddress] = accumulator
            }
          }
        }
      }
    }
  }
}
