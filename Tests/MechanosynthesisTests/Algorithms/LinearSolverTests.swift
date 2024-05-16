import XCTest
import Mechanosynthesis
import Numerics

// ========================================================================== //
// Performance of different algorithms
// ========================================================================== //
//
// Specification of multigrid V-cycle:
//   A-B-C-D-E
//
//   A GSRB (1x)
//   -> B GSRB (2x)
//   ---> C GSRB (4x)
//   -> D GSRB (2x)
//   E GSRB (1x)
//   update solution
//
//   One V-cycle counts as one iteration. It is effectively the compute cost
//   of two Gauss-Seidel iterations.
//
// ========================================================================== //
// Raw data
// ========================================================================== //
//
// h = 0.25, gridSize = 8, cellCount = 512
//                       0 iters  ||r|| = 394.27557
// Gauss-Seidel         30 iters  ||r|| = 6.5650797      0.002 seconds
// Conjugate Gradient   30 iters  ||r|| = 0.00030688613  0.003 seconds
// Multigrid 1-1-1-1-1  15 iters  ||r|| = 0.017259505    0.004 seconds
// Preconditioned CG    15 iters  ||r|| = 0.00023875975  0.003 seconds
//
// h = 0.125, gridSize = 16, cellCount = 4096
//                       0 iters  ||r|| = 3091.9424
// Gauss-Seidel         30 iters  ||r|| = 277.43747     0.016 seconds
// Conjugate Gradient   30 iters  ||r|| = 0.09551496    0.017 seconds
// Multigrid 1-1-1-1-1  15 iters  ||r|| = 0.075272486   0.033 seconds
// Preconditioned CG    15 iters  ||r|| = 0.0032440922  0.024 seconds
//
// h = 0.0625, gridSize = 32, cellCount = 32,768
//                       0 iters  ||r|| = 24494.229
// Gauss-Seidel         60 iters  ||r|| = 1308.8044    0.250 seconds
// Conjugate Gradient   60 iters  ||r|| = 0.49065304   0.258 seconds
// Multigrid 1-2-2-2-1  30 iters  ||r|| = 0.035719506  0.542 seconds
// Preconditioned CG    30 iters  ||r|| = 0.048568394  0.364 seconds
//
// h = 0.0313, gridSize = 64, cellCount = 262,144
//                       0 iters  ||r|| = 195015.61
// Gauss-Seidel         99 iters  ||r|| = 5887.104   3.311 seconds
// Conjugate Gradient   99 iters  ||r|| = 53.441914  3.375 seconds
// Multigrid 1-2-4-2-1  60 iters  ||r|| = 0.3711439  8.855 seconds
// Preconditioned CG    60 iters  ||r|| = 0.7311449  5.874 seconds
//
// h = 0.0156, gridSize = 128, cellCount = 2,097,152
// Multigrid 1-2-2-1-2-2-1  60 iters  ||r|| = 294.65112  70.652 seconds
// Multigrid 1-2-2-2-2-2-1  60 iters  ||r|| = 38.89426   71.247 seconds
// Multigrid 1-2-2-4-2-2-1  60 iters  ||r|| = 3.4102564  72.914 seconds
// Preconditioned CG        60 iters  ||r|| = 1209.9086  46.300 seconds
// Preconditioned CG        99 iters  ||r|| = 11.659912  75.554 seconds
//
// ========================================================================== //
// Conclusions
// ========================================================================== //
//
// Ranked in order of ease of implementation:
// 1) Jacobi
// 2) Gauss-Seidel
// 3) Conjugate gradient
// 4) Preconditioned conjugate gradient
// 5) Multigrid
//
// Preconditioned CG seems like the best tradeoff between complexity and speed.
// It converges consistently in every situation. Multigrid requires careful
// tuning of the tree depth and often fails to converge with the wrong V-cycle
// scheme. However, it has the potential to be more efficient, especially with
// BMR FAS-FMG.
//
// I'm also unsure how adaptive mesh refinement will affect the performance of
// these algorithms. The path length to jump between levels would increase
// significantly. Multigrid would coalesce the overhead of interpolation
// and coarsening operations. However, the CG preconditioner could be modified
// with higher / anisotropic sample count at the nuclear singularity.
//
// The idea of multigrid F-cycles could be transferred over to CG. In addition,
// I'd like to consider a hybrid between multigrid and preconditioned CG.
// Multigrid performs robustly whenever there's only 2 levels. Perhaps one
// could start a PCG loop within the lower resolution level, which lasts for
// around 8 iterations.
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
    for _ in 0..<60 {
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
    var Kr = applyLaplacianPreconditioner(r)
    var rKr = Self.dot(r, Kr)
    var p = Kr
    
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
      
      let Ap = Self.applyLaplacianLinearPart(p)
      
      let a = rKr / Self.dot(p, Ap)
      let xNew = Self.shift(x, scale: a, correction: p)
      let rNew = Self.shift(r, scale: -a, correction: Ap)
      let KrNew = applyLaplacianPreconditioner(rNew)
      let rKrNew = Self.dot(rNew, KrNew)
      
      let b = rKrNew / rKr
      let pNew = Self.shift(KrNew, scale: b, correction: p)
      
      x = xNew
      r = rNew
      Kr = KrNew
      rKr = rKrNew
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
  //
  // a four-color scheme would work with Mehrstellen, provided we process the
  // multigrid one level at a time
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
    
    print()
    print("Gauss-Seidel")
    for iterationID in 0..<30 {
      do {
        let L1x = Self.applyLaplacianLinearPart(x)
        let r = Self.shift(b, scale: -1, correction: L1x)
        let r2 = Self.dot(r, r)
        let normres = r2.squareRoot()
        print("||r|| = \(normres)")
      }
      
      executeSweep(red: true, black: false)
      executeSweep(red: false, black: true)
    }
  }
  
  func testMultigridMethod() {
    var b = Self.createScaledChargeDensity()
    var x = [Float](repeating: .zero, count: Self.cellCount)
    
    // One V-cycle should be treated as two SD or CG iterations.
    print()
    print("Multigrid")
    for iterationID in 0..<15 {
      do {
        let L1x = Self.applyLaplacianLinearPart(x)
        let L2x = Self.applyLaplacianBoundary()
        let Ax = Self.shift(L1x, scale: 1, correction: L2x)
        let r = Self.shift(b, scale: -1, correction: Ax)
        let r2 = Self.dot(r, r)
        let normres = r2.squareRoot()
        print("||r|| = \(normres)")
      }
      
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
      GSRB_LEVEL(e: &eFine, r: rFine, coarseness: 1, iteration: 0)
      GSRB_LEVEL(e: &eFine, r: rFine, coarseness: 1, iteration: 1)
      multigridCoarseLevel(eFine: &eFine, rFine: &rFine, fineLevelCoarseness: 1)
      
      // Update the solution.
      x = Self.shift(x, scale: 1, correction: eFine)
    }
    
    func multigridCoarseLevel(
      eFine: inout [Float], rFine: inout [Float], fineLevelCoarseness: Int
    ) {
      // Determine the array lengths.
      var fineArrayLength = Self.cellCount
      fineArrayLength /= fineLevelCoarseness
      fineArrayLength /= fineLevelCoarseness
      fineArrayLength /= fineLevelCoarseness
      
      // Restrict from fine to coarse.
      var rFineCorrected = rFine
      var rCoarse: [Float] = []
      correctResidual(
        e: eFine, 
        r: &rFineCorrected,
        coarseness: fineLevelCoarseness)
      shiftResolution(
        fineGrid: &rFineCorrected,
        coarseGrid: &rCoarse,
        fineLevelCoarseness: fineLevelCoarseness, 
        shiftingUp: true)
      
      // Smoothing iterations on the coarse level.
      let coarseLevelCoarseness = 2 * fineLevelCoarseness
      var eCoarse = [Float](repeating: .zero, count: fineArrayLength / 8)
      GSRB_LEVEL(
        e: &eCoarse, 
        r: rCoarse,
        coarseness: coarseLevelCoarseness,
        iteration: 0)
      GSRB_LEVEL(
        e: &eCoarse,
        r: rCoarse,
        coarseness: coarseLevelCoarseness,
        iteration: 1)
      
      // Shift to a higher level.
      if coarseLevelCoarseness < 4 {
        multigridCoarseLevel(
          eFine: &eCoarse,
          rFine: &rCoarse,
          fineLevelCoarseness: coarseLevelCoarseness)
      }
      
      // Prolong from coarse to fine.
      shiftResolution(
        fineGrid: &eFine, 
        coarseGrid: &eCoarse,
        fineLevelCoarseness: fineLevelCoarseness, 
        shiftingUp: false)
      correctResidual(
        e: eFine,
        r: &rFine,
        coarseness: fineLevelCoarseness)
      
      // Smoothing iterations on the fine level.
      var δeFine = [Float](repeating: .zero, count: fineArrayLength)
      GSRB_LEVEL(
        e: &δeFine, 
        r: rFine,
        coarseness: fineLevelCoarseness,
        iteration: 0)
      GSRB_LEVEL(
        e: &δeFine,
        r: rFine,
        coarseness: fineLevelCoarseness, 
        iteration: 1)
      for cellID in 0..<fineArrayLength {
        eFine[cellID] += δeFine[cellID]
      }
    }
    
    // Gauss-Seidel red-black
    func GSRB_LEVEL(
      e: inout [Float], r: [Float], coarseness: Int, iteration: Int
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
            guard (iteration & 1) == (parity & 1) else {
              continue
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
              if iteration == 0 {
                let neighborValue = r[neighborAddress]
                let λ = h * h / 6
                faceAccumulator += 1 / (h * h) * (-λ * neighborValue)
              } else {
                let neighborValue = e[neighborAddress]
                faceAccumulator += 1 / (h * h) * neighborValue
              }
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
  
  // A hybrid between multigrid and conjugate gradient.
  // - CG nested inside of MG.
  // - MG only spans 2 levels.
  // - Using the 33-point convolution to precondition CG.
  func testHybridMethod() throws {
    // TODO: Optimize the multigrid method before creating the hybrid kernel.
    // It may have been implemented incorrectly, explaining some of the
    // instability.
  }
}
