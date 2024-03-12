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
  
  /// The block size for intermediate band reduction.
  ///
  /// If not specified, the optimal block size for the hardware backend is
  /// chosen.
  public var blockSize: Int?
  
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
  var matrix: [Float]
  
  /// The number of unknowns, `n`.
  var problemSize: Int
  
  /// The panel factorization size, `nb`.
  var blockSize: Int
  
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
    self.matrix = matrix
    self.problemSize = problemSize
    self.blockSize = -1
    
    guard matrix.count == problemSize * problemSize else {
      fatalError("""
        Invalid matrix size: expected '\(problemSize * problemSize)' \
        but got '\(matrix.count)'.
        """)
    }
    guard problemSize > 0 else {
      fatalError("Cannot solve a problem with less than one unknown.")
    }
    guard problemSize > 1 else {
      // This cannot be solved with standard tridiagonalization techniques.
      eigenvalues = [matrix[0]]
      eigenvectors = [1]
      return
    }
    
    createBlockSize(descriptor: descriptor)
    solveEigenproblem()
  }
  
  mutating func createBlockSize(descriptor: DiagonalizationDescriptor) {
    if let blockSize = descriptor.blockSize {
      guard blockSize > 0 else {
        fatalError("Block size must be at least one.")
      }
      self.blockSize = blockSize
    } else {
      self.blockSize = 4
    }
    blockSize = min(blockSize, problemSize - 1)
  }
  
  mutating func solveEigenproblem() {
    let bandFormReflectors = reduceToBandForm()
    let bulgeChasingReflectors = chaseBulges()
    solveTridiagonalEigenproblem()
    backTransform(
      bandFormReflectors: bandFormReflectors,
      bulgeChasingReflectors: bulgeChasingReflectors)
    
    // Deallocate the matrix after it's finished.
    matrix = []
  }
}
