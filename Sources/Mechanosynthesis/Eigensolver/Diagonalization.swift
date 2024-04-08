//
//  Diagonalize.swift
//
//
//  Created by Philip Turner on 3/11/24.
//

import QuartzCore

/// A configuration for an eigendecomposition.
public struct DiagonalizationDescriptor {
  /// A symmetric matrix of FP32 numbers.
  ///
  /// The size must equal `n` squared.
  public var matrix: [Float]?
  
  /// The number of unknowns to solve for.
  public var problemSize: Int?
  
  /// The block size for panel factorization.
  ///
  /// If not specified, the optimal block size for the hardware backend is
  /// chosen.
  public var blockSize: Int?
  
  /// The block size for recursive panel factorization.
  ///
  /// If not specified, the optimal block size for the hardware backend is
  /// chosen.
  public var smallBlockSize: Int?
  
  public init() {
    
  }
}

/// The decomposition of a matrix into its principal components.
public struct Diagonalization {
  /// The eigenvalues, in ascending order.
  ///
  /// The number of array elements equals `n`. In the case of degenerate
  /// eigenvalue clusters, an eigenvalue may appear multiple times.
  public var eigenvalues: [Float] = []
  
  /// A column-major matrix of eigenvectors.
  ///
  /// Each eigenvector corresponds to the eigenvalue at the same position.
  public var eigenvectors: [Float] = []
  
  /// The current state of the matrix, which is progressing toward tridiagonal
  /// form.
  var matrix: [Float] = []
  var matrixPointer: UnsafeMutablePointer<Float>?
  
  /// The number of unknowns, `n`.
  var problemSize: Int
  
  /// The panel factorization size, `nb`.
  var blockSize: Int
  
  /// The recursive panel factorization size, `sb`.
  var smallBlockSize: Int
  
  /// Form the eigendecomposition of a matrix.
  ///
  /// This is a blocking function. It immediately solves the eigenproblem,
  /// which may incur a large amount of latency. It may also perform a
  /// blocking access to GPU hardware.
  public init(descriptor: DiagonalizationDescriptor) {
    
    guard let matrix = descriptor.matrix,
          let problemSize = descriptor.problemSize else {
      fatalError("Invalid descriptor.")
    }
    guard matrix.count == problemSize * problemSize else {
      fatalError("""
        Invalid matrix size: expected '\(problemSize * problemSize)' \
        but got '\(matrix.count)'.
        """)
    }
    guard problemSize > 0 else {
      fatalError("Cannot solve a problem with less than one unknown.")
    }
    self.problemSize = problemSize
    self.blockSize = -1
    self.smallBlockSize = -1
    
    guard problemSize > 1 else {
      // This cannot be solved with standard tridiagonalization techniques.
      eigenvalues = [matrix[0]]
      eigenvectors = [1]
      return
    }
    
    
    createBlockSize(descriptor: descriptor)
    
    solveEigenproblem(descriptor: descriptor)
   
  }
  
  mutating func createBlockSize(descriptor: DiagonalizationDescriptor) {
    // Switch over the permutations of available pieces of information.
    switch (descriptor.blockSize, descriptor.smallBlockSize) {
    case (.some(let blockSize), .some(let smallBlockSize)):
      self.blockSize = blockSize
      self.smallBlockSize = smallBlockSize
    case (.some(let blockSize), .none):
      self.blockSize = blockSize
      self.smallBlockSize = max(1, min((blockSize + 3) / 4, blockSize))
    case (.none, .some(let smallBlockSize)):
      fatalError("""
        The small block size was specified (\(smallBlockSize)), \
        but the large block size was not.
        """)
    case (.none, .none):
      // The current heuristic is a placeholder for a proper one.
      blockSize = 4
      smallBlockSize = 2
    }
    
    guard blockSize <= problemSize else {
      fatalError("Block size cannot exceed problem size.")
    }
    guard blockSize > 0, smallBlockSize > 0 else {
      fatalError("Block size must be at least one.")
    }
    guard blockSize % smallBlockSize == 0 else {
      fatalError("Block size must be divisible by small block size.")
    }
  }
  
  mutating func createMatrix(descriptor: DiagonalizationDescriptor) {
    matrix = Array(repeating: .zero, count: problemSize * problemSize)
    memcpy(&matrix, descriptor.matrix!, problemSize * problemSize * 4)
    matrixPointer = matrix.withContiguousMutableStorageIfAvailable {
      $0.baseAddress!
    }!
    guard UInt(bitPattern: matrixPointer) != 0 else {
      fatalError("Matrix pointer was null.")
    }
  }
  
  // Deallocate the matrix after it's finished.
  mutating func destroyMatrix() {
    matrix = []
    matrixPointer = nil
  }
  
  mutating func solveEigenproblem(descriptor: DiagonalizationDescriptor) {
    var bandReflectors: [Float] = []
    var bulgeReflectors: [Float] = []
    
    // Reduce the bandwidth from 'problemSize' to 'blockSize'.
    let checkpoint0 = CACurrentMediaTime()
    createMatrix(descriptor: descriptor)
    bandReflectors = reduceToBandForm()
    
    // If the matrix is already in tridiagonal form, there is no work to do.
    let checkpoint1 = CACurrentMediaTime()
    if blockSize > 1 {
      bulgeReflectors = chaseBulges()
    }
    
    // Acquire the eigenvectors (using the LAPACK divide-and-conquer until all
    // other bottlenecks are suppressed).
    let checkpoint2 = CACurrentMediaTime()
    solveTridiagonalEigenproblem()
    
    // If the matrix was already in tridiagonal form, there is no work to do.
    let checkpoint3 = CACurrentMediaTime()
    if blockSize > 1 {
      backTransform(bulgeReflectors: bulgeReflectors)
    }
    
    // Expand the bandwidth from 'blockSize' to 'problemSize'.
    let checkpoint4 = CACurrentMediaTime()
    backTransform(bandReflectors: bandReflectors)
    destroyMatrix()
    let checkpoint5 = CACurrentMediaTime()
    
    // Report diagonistics for debugging performance.
    if problemSize >= 100 {
      let time01 = 1e6 * (checkpoint1 - checkpoint0)
      let time12 = 1e6 * (checkpoint2 - checkpoint1)
      let time23 = 1e6 * (checkpoint3 - checkpoint2)
      let time34 = 1e6 * (checkpoint4 - checkpoint3)
      let time45 = 1e6 * (checkpoint5 - checkpoint4)
      let times = SIMD8(time01, time12, time23, time34, time45, 0, 0, 0)
      
      let total = times.sum()
      let percents = times / total * 100
      func format(_ number: Double) -> String {
        String(format: "%.1f", number)
      }
      func printPart(index: Int, label: String) {
        print("[\(format(percents[index]))%", terminator: ", ")
        print(format(times[index]), terminator: " μs] - ")
        print(label)
      }
      
      print()
      print("[n = \(problemSize)] Performance Breakdown:")
      printPart(index: 0, label: "Reduction to band form")
      printPart(index: 1, label: "Bulge chasing")
      printPart(index: 2, label: "Divide and conquer")
      printPart(index: 3, label: "Back transformation (2nd stage)")
      printPart(index: 4, label: "Back transformation (1st stage)")
    }
  }
}
