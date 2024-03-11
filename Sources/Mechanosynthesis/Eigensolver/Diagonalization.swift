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
  public var n: Int?
  
  /// The block size for intermediate band reduction.
  ///
  /// If not specified, the optimal block size for the hardware backend is
  /// chosen.
  public var nb: Int?
  
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
          let n = descriptor.n else {
      fatalError("Invalid descriptor.")
    }
    self.matrix = matrix
    self.problemSize = n
    self.blockSize = -1
    
    guard matrix.count == n * n else {
      fatalError(
        "Invalid matrix size: expected '\(n * n)' but got '\(matrix.count)'.")
    }
    guard n > 0 else {
      fatalError("Cannot solve a problem with less than one unknown.")
    }
    guard n > 1 else {
      // This cannot be solved with standard tridiagonalization techniques.
      eigenvalues = [matrix[0]]
      eigenvectors = [1]
      return
    }
    
    createBlockSize(descriptor: descriptor)
    solveEigenproblem()
  }
  
  mutating func createBlockSize(descriptor: DiagonalizationDescriptor) {
    if let nb = descriptor.nb {
      guard nb > 0 else {
        fatalError("Block size must be at least one.")
      }
      blockSize = nb
    } else {
      blockSize = 4
    }
    blockSize = min(blockSize, problemSize - 1)
  }
  
  mutating func solveEigenproblem() {
    (eigenvalues, eigenvectors) = diagonalize(
      matrix: matrix, n: problemSize, nb: blockSize)
    matrix = []
  }
}

// Decomposes a matrix into its principal components.
//
// Arguments:
// - matrix: symmetric matrix of FP32 numbers.
// - n: number of unknowns to solve for.
// - nb: block size for intermediate band reduction.
//
// Returns:
// - eigenvalues: n-element array of eigenvalues, in ascending order.
// - eigenvectors: column-major matrix of the associated eigenvectors.
func diagonalize(matrix: [Float], n: Int, nb: Int) -> (
  eigenvalues: [Float], eigenvectors: [Float]
) {
  // Allocate main memory allocations.
  var currentMatrixA = matrix
  var currentReflectors = [Float](repeating: 0, count: n * n)
  
  // Reduce the matrix to band form, and collect up the reflectors.
  var blockStart: Int = 0
  while blockStart < n - nb {
    // Adjust the loop end, to account for the factorization band offset.
    let blockEnd = min(blockStart + nb, n - nb)
    defer { blockStart += nb }
    
    // Load to panel into the cache, isolating mutations to the matrix A.
    var panel = [Float](repeating: 0, count: nb * n)
    for rowID in blockStart..<blockEnd {
      for columnID in 0..<n {
        let matrixAddress = rowID * n + columnID
        let panelAddress = (rowID - blockStart) * n + columnID
        panel[panelAddress] = currentMatrixA[matrixAddress]
      }
    }
    
    // Allocate cache memory for the reflectors.
    var panelReflectors = [Float](repeating: 0, count: nb * n)
    
    // Generate the reflectors.
    for reflectorID in blockStart..<blockEnd {
      // Factor starting at an offset from the diagonal.
      let bandOffset = reflectorID + nb
      
      // Load the row into the cache.
      var vector = [Float](repeating: 0, count: n)
      for elementID in 0..<n {
        let address = (reflectorID - blockStart) * n + elementID
        vector[elementID] = panel[address]
      }
      
      // Apply preceding reflectors (from this panel) to the column.
      for previousReflectorID in blockStart..<reflectorID {
        // Load the reflector into the cache.
        var reflector = [Float](repeating: 0, count: n)
        for elementID in 0..<n {
          let address = (previousReflectorID - blockStart) * n + elementID
          reflector[elementID] = panelReflectors[address]
        }
        
        // Apply the reflector.
        var dotProduct: Float = .zero
        for elementID in 0..<n {
          dotProduct += reflector[elementID] * vector[elementID]
        }
        for elementID in 0..<n {
          vector[elementID] -= reflector[elementID] * dotProduct
        }
      }
      
      // Zero out the elements above the band offset.
      for elementID in 0..<bandOffset {
        vector[elementID] = 0
      }
      
      // Take the norm of the vector.
      var norm: Float = .zero
      for elementID in 0..<n {
        norm += vector[elementID] * vector[elementID]
      }
      norm.formSquareRoot()
      
      // Modify the vector, turning it into a reflector.
      let oldSubdiagonal = vector[bandOffset]
      let newSubdiagonal = norm * Float((oldSubdiagonal >= 0) ? -1 : 1)
      let tau = (newSubdiagonal - oldSubdiagonal) / newSubdiagonal
      for elementID in 0..<n {
        var element = vector[elementID]
        if elementID == bandOffset {
          element = 1
        } else {
          element /= oldSubdiagonal - newSubdiagonal
        }
        element *= tau.squareRoot()
        vector[elementID] = element
      }
      
      // Store the reflector to the cache.
      for elementID in 0..<n {
        let address = (reflectorID - blockStart) * n + elementID
        panelReflectors[address] = vector[elementID]
      }
    }
    
    // Apply the reflectors to the matrix, from both sides.
    for directionID in 0..<2 {
      for vectorID in 0..<n {
        var vector = [Float](repeating: 0, count: n)
        if directionID == 0 {
          // Load the row into the cache.
          for elementID in 0..<n {
            let address = vectorID * n + elementID
            vector[elementID] = currentMatrixA[address]
          }
        } else {
          // Load the column into the cache.
          for elementID in 0..<n {
            let address = elementID * n + vectorID
            vector[elementID] = currentMatrixA[address]
          }
        }
        
        for reflectorID in blockStart..<blockEnd {
          // Load the reflector into the cache.
          var reflector = [Float](repeating: 0, count: n)
          for elementID in 0..<n {
            let address = (reflectorID - blockStart) * n + elementID
            reflector[elementID] = panelReflectors[address]
          }
          
          // Apply the reflector.
          var dotProduct: Float = .zero
          for elementID in 0..<n {
            dotProduct += reflector[elementID] * vector[elementID]
          }
          for elementID in 0..<n {
            vector[elementID] -= reflector[elementID] * dotProduct
          }
        }
        
        if directionID == 0 {
          // Store the row to main memory.
          for elementID in 0..<n {
            let address = vectorID * n + elementID
            currentMatrixA[address] = vector[elementID]
          }
        } else {
          // Store the column to main memory.
          for elementID in 0..<n {
            let address = elementID * n + vectorID
            currentMatrixA[address] = vector[elementID]
          }
        }
      }
    }
    
    // Store the reflectors to main memory.
    for reflectorID in blockStart..<blockEnd {
      for elementID in 0..<n {
        let cacheAddress = (reflectorID - blockStart) * n + elementID
        let memoryAddress = reflectorID * n + elementID
        currentReflectors[memoryAddress] = panelReflectors[cacheAddress]
      }
    }
  }
  
  // MARK: - Bulge Chasing
  
  var bulgeReflectors: [[Float]] = []
  
  // Apply the Householder reflector to the entire matrix. This isn't very
  // efficient, but it's correct, which we want for initial debugging.
  func applyBulgeChase(
    sweepID: Int,
    vectorID: Int,
    startElementID: Int,
    endElementID: Int
  ) {
    // Load the row into the cache.
    var vector = [Float](repeating: 0, count: n)
    for elementID in startElementID..<endElementID {
      let address = vectorID * n + elementID
      vector[elementID] = currentMatrixA[address]
    }
    
    // Take the norm of the vector.
    var norm: Float = .zero
    for elementID in startElementID..<endElementID {
      norm += vector[elementID] * vector[elementID]
    }
    norm.formSquareRoot()
    
    // Modify the vector, turning it into a reflector.
    let oldSubdiagonal = vector[startElementID]
    let newSubdiagonal = norm * Float((oldSubdiagonal >= 0) ? -1 : 1)
    let tau = (newSubdiagonal - oldSubdiagonal) / newSubdiagonal
    for elementID in startElementID..<endElementID {
      var element = vector[elementID]
      if elementID == startElementID {
        element = 1
      } else {
        element /= oldSubdiagonal - newSubdiagonal
      }
      element *= tau.squareRoot()
      vector[elementID] = element
    }
    
    // Apply the reflector to the matrix, from both sides.
    let reflector = vector
    for directionID in 0..<2 {
      for vectorID in 0..<n {
        var vector = [Float](repeating: 0, count: n)
        if directionID == 0 {
          // Load the row into the cache.
          for elementID in 0..<n {
            let address = vectorID * n + elementID
            vector[elementID] = currentMatrixA[address]
          }
        } else {
          // Load the column into the cache.
          for elementID in 0..<n {
            let address = elementID * n + vectorID
            vector[elementID] = currentMatrixA[address]
          }
        }
        
        // Apply the reflector.
        var dotProduct: Float = .zero
        for elementID in 0..<n {
          dotProduct += reflector[elementID] * vector[elementID]
        }
        for elementID in 0..<n {
          vector[elementID] -= reflector[elementID] * dotProduct
        }
        
        if directionID == 0 {
          // Store the row to main memory.
          for elementID in 0..<n {
            let address = vectorID * n + elementID
            currentMatrixA[address] = vector[elementID]
          }
        } else {
          // Store the column to main memory.
          for elementID in 0..<n {
            let address = elementID * n + vectorID
            currentMatrixA[address] = vector[elementID]
          }
        }
      }
    }
    
    // Store the reflector to main memory.
    bulgeReflectors.append(reflector)
  }
  
  if nb > 1 {
    for sweepID in 0..<max(0, n - 2) {
      let startVectorID = sweepID + 1
      var endVectorID = sweepID + nb + 1
      endVectorID = min(endVectorID, n)
      guard endVectorID - startVectorID > 1 else {
        fatalError("Generated empty Householder transform.")
      }
      
      applyBulgeChase(
        sweepID: sweepID, vectorID: sweepID,
        startElementID: startVectorID, endElementID: endVectorID)
      
      var operationID = 1
      while true {
        let startColumnID = (sweepID - nb + 1) + operationID * nb
        let startVectorID = (sweepID + 1) + operationID * nb
        var endVectorID = (sweepID + nb + 1) + operationID * nb
        endVectorID = min(endVectorID, n)
        
        if endVectorID - startVectorID > 1 {
          applyBulgeChase(
            sweepID: sweepID, vectorID: startColumnID,
            startElementID: startVectorID, endElementID: endVectorID)
        } else {
          break
        }
        operationID += 1
      }
    }
  }
  
  // MARK: - Validation Testing
  
  var (eigenvalues, eigenvectors) = divideAndConquer(
    matrix: currentMatrixA, n: n)
  
  // Back-transform the eigenvectors.
  for vectorID in 0..<n {
    // Load the vector into the cache.
    var vector = [Float](repeating: 0, count: n)
    for elementID in 0..<n {
      let address = vectorID * n + elementID
      vector[elementID] = eigenvectors[address]
    }
    
    for reflectorID in bulgeReflectors.indices.reversed() {
      // Load the reflector into the cache.
      let reflector = bulgeReflectors[reflectorID]
      
      // Apply the reflector.
      var dotProduct: Float = .zero
      for elementID in 0..<n {
        dotProduct += reflector[elementID] * vector[elementID]
      }
      for elementID in 0..<n {
        vector[elementID] -= reflector[elementID] * dotProduct
      }
    }
    
    for reflectorID in (0..<n).reversed() {
      // Load the reflector into the cache.
      var reflector = [Float](repeating: 0, count: n)
      for elementID in 0..<n {
        let address = reflectorID * n + elementID
        reflector[elementID] = currentReflectors[address]
      }
      
      // Apply the reflector.
      var dotProduct: Float = .zero
      for elementID in 0..<n {
        dotProduct += reflector[elementID] * vector[elementID]
      }
      for elementID in 0..<n {
        vector[elementID] -= reflector[elementID] * dotProduct
      }
    }
    
    // Store the vector to main memory.
    for elementID in 0..<n {
      let address = vectorID * n + elementID
      eigenvectors[address] = vector[elementID]
    }
  }
  
  return (eigenvalues, eigenvectors)
}
