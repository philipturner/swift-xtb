//
//  Diagonalize.swift
//
//
//  Created by Philip Turner on 3/11/24.
//

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
      if problemSize >= 8 {
        blockSize = min(32, problemSize & (~7))
        smallBlockSize = 8
      } else {
        blockSize = problemSize
        smallBlockSize = 1
      }
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
    // Allocate a new array.
    let elementCount = problemSize * problemSize
    matrix = Array(repeating: .zero, count: elementCount)
    
    // Copy the elements over, avoiding unsafe behavior from CoW.
    for elementID in 0..<elementCount {
      let descriptorElement = descriptor.matrix.unsafelyUnwrapped[elementID]
      matrix[elementID] = descriptorElement
    }
    
    // Retrieve the pointer and use it for the remaining computations.
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
    createMatrix(descriptor: descriptor)
    bandReflectors = reduceToBandForm()

    // If the matrix is already in tridiagonal form, there is no work to do.
    if blockSize > 1 {
      bulgeReflectors = chaseBulges()
    }

    // Acquire the eigenvectors (using the LAPACK divide-and-conquer until all
    // other bottlenecks are suppressed).
    solveTridiagonalEigenproblem()

    // If the matrix was already in tridiagonal form, there is no work to do.
    if blockSize > 1 {
      backTransform(bulgeReflectors: bulgeReflectors)
    }

    // Expand the bandwidth from 'blockSize' to 'problemSize'.
    backTransform(bandReflectors: bandReflectors)
    destroyMatrix()
  }
}
