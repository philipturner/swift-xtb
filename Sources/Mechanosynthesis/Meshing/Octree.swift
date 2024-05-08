//
//  Octree.swift
//
//
//  Created by Philip Turner on 2/18/24.
//

/// A configuration for an octree.
public struct OctreeDescriptor {
  /// Required. The power-2 size of the coarsest level.
  public var sizeExponent: Int?
  
  public init() {
    
  }
}

/// The data needed for locating a node in 3D space and memory.
public struct OctreeNode: Equatable {
  /// A four-element vector that compactly stores the center and spacing.
  public var centerAndSpacing: SIMD4<Float>
  
  /// The average 3D position of the cell's contents.
  @_transparent public var center: SIMD3<Float> {
    unsafeBitCast(centerAndSpacing, to: SIMD3<Float>.self)
  }
  
  /// The cube root of the volume.
  @_transparent public var spacing: Float {
    centerAndSpacing.w
  }
  
  /// The array index of the parent. If there is none, this is 'UInt32.max'.
  public var parentIndex: UInt32
  
  /// The array index where the child nodes start. If there are none, this is
  /// 'UInt32.max'.
  public var branchesIndex: UInt32
  
  /// If any children have grand-children, they are marked in this mask. Such
  /// children are stored in compacted order. If the child is a leaf node, it
  /// isn't marked.
  public var branchesMask: SIMD8<UInt8>
  
  @_transparent
  public init(
    centerAndSpacing: SIMD4<Float>,
    parentIndex: UInt32,
    branchesIndex: UInt32,
    branchesMask: SIMD8<UInt8>
  ) {
    self.centerAndSpacing = centerAndSpacing
    self.parentIndex = parentIndex
    self.branchesIndex = branchesIndex
    self.branchesMask = branchesMask
  }
}

/// An octree data structure designed for efficient traversal.
public struct Octree {
  /// The cells from every hierarchy level, in Morton order.
  public var nodes: [OctreeNode] = []
  
  public init(descriptor: OctreeDescriptor) {
    guard let sizeExponent = descriptor.sizeExponent else {
      fatalError("Descriptor was invalid.")
    }
    
    let origin: SIMD3<Float> = .zero
    let size = Float(sign: .plus, exponent: sizeExponent, significand: 1)
    let node = OctreeNode(
      centerAndSpacing: SIMD4(origin, size),
      parentIndex: UInt32.max,
      branchesIndex: UInt32.max,
      branchesMask: SIMD8(repeating: UInt8.max))
    nodes = [node]
  }
  
  /// Efficient function to expand/contract several nodes at once.
  /// - Parameter expanded: The nodes that should acquire children. Each array
  ///   element is the parent address and the mask of children that should be
  ///   expanded.
  /// - Parameter contracted: The parents who should no longer have children.
  ///   The children cannot have grandchildren.
  /// - returns: The old nodes' positions in the new list. If no such position
  ///            exists, the array element is `UInt32.max`.
  @discardableResult
  public mutating func resizeNodes(
    expanded insertionMarks: [SIMD8<UInt8>],
    contracted removalMarks: [Bool]
  ) -> [UInt32] {
    @_transparent
    func isMax(_ vector: SIMD8<UInt8>) -> Bool {
      unsafeBitCast(vector, to: UInt64.self) == .max
    }
    
    // Create an array of arrays for the hierarchy levels.
    var levels: [[OctreeNode]] = []
    var mappedIndices: [[UInt32]] = []
    levels.append([])
    mappedIndices.append([])
    
    func traversalFirstPass(nodeID: Int, levelID: Int) {
      var node = nodes[nodeID]
      if removalMarks[nodeID] {
        fatalError("Traversed to a node that should be removed.")
      }
      if levelID != 0 {
        node.parentIndex = UInt32(levels[levelID - 1].count)
      } else {
        node.parentIndex = UInt32.max
      }
      if levels.count <= levelID + 1 {
        levels.append([])
        mappedIndices.append([])
      }
      
      var previousBranchesIndex = Int(node.branchesIndex)
      let nextBranchesIndex = UInt32(levels[levelID + 1].count)
      
      // Lookup table for child nodes.
      var x = SIMD8<Float>(0, 1, 0, 1, 0, 1, 0, 1) * 0.5 - 0.25
      var y = SIMD8<Float>(0, 0, 1, 1, 0, 0, 1, 1) * 0.5 - 0.25
      var z = SIMD8<Float>(0, 0, 0, 0, 1, 1, 1, 1) * 0.5 - 0.25
      let selfIndex = UInt32(levels[levelID].count)
      x = x * node.spacing + node.center.x
      y = y * node.spacing + node.center.y
      z = z * node.spacing + node.center.z
      
      var nodeBranchesMask = unsafeBitCast(node.branchesMask, to: UInt64.self)
      let insertion = unsafeBitCast(insertionMarks[nodeID], to: UInt64.self)
      
      for branchID in 0..<8 {
        let branchBits = UInt64(255) << UInt64(branchID * 8)
        
        if nodeBranchesMask & branchBits != branchBits {
          if removalMarks[previousBranchesIndex] {
            nodeBranchesMask |= branchBits
            mappedIndices[levelID + 1].append(UInt32.max)
          } else {
            traversalFirstPass(
              nodeID: previousBranchesIndex, levelID: levelID + 1)
          }
          previousBranchesIndex += 1
        } else if insertion & branchBits != 0 {
          nodeBranchesMask &= ~branchBits
          
          // Construct child position using the lookup table.
          let xyz = SIMD3(x[branchID], y[branchID], z[branchID])
          
          // Insert a node for the child.
          let node = OctreeNode(
            centerAndSpacing: SIMD4(xyz, node.spacing / 2),
            parentIndex: UInt32(selfIndex),
            branchesIndex: UInt32.max,
            branchesMask: SIMD8(repeating: UInt8.max))
          levels[levelID + 1].append(node)
        }
      }
      
      node.branchesMask = unsafeBitCast(nodeBranchesMask, to: SIMD8<_>.self)
      if isMax(node.branchesMask) {
        node.branchesIndex = UInt32.max
      } else {
        node.branchesIndex = UInt32(nextBranchesIndex)
        
        func prefixSum(_ x: SIMD8<UInt8>) -> SIMD8<UInt8> {
          // 00000000
          // 11111111
          var output: UInt64 = .zero
          let forwardSum0 = unsafeBitCast(x, to: UInt64.self)
          
          // 01010101
          // 02020202
          let shifted1 = (forwardSum0 & 0x00FF00FF00FF00FF) << 8
          let forwardSum1 = (forwardSum0 &+ shifted1) & 0xFF00FF00FF00FF00
          output = shifted1
          
          // 01230123
          // 00040004
          let shifted2 = (forwardSum1 & 0x0000FF000000FF00) << 16
          let forwardSum2 = (forwardSum1 &+ shifted2) & 0xFF000000FF000000
          output += (forwardSum1 & 0x0000FF000000FF00) * (0x0101 << 8)
          
          // 01234567
          output += forwardSum2 * (0x01010101) << 8
          return unsafeBitCast(output, to: SIMD8<UInt8>.self)
        }
        
        var x = SIMD8<UInt8>(repeating: 0)
        x.replace(with: 1, where: node.branchesMask .!= 255)
        x = prefixSum(x)
        node.branchesMask.replace(with: x, where: node.branchesMask .!= 255)
      }
      mappedIndices[levelID].append(UInt32(selfIndex))
      levels[levelID].append(node)
    }
    
    // Invoke the recursive function from the top level.
    traversalFirstPass(nodeID: 0, levelID: 0)
    
    // Bookkeeping operations on the indices.
    let mappedIndexCount = mappedIndices.reduce(0) { $0 + $1.count }
    guard mappedIndexCount == nodes.count else {
      fatalError("Unexpected behavior when mapping old nodes to new nodes.")
    }
    var levelPrefixSums: [UInt32] = []
    var levelSum: Int = .zero
    for level in levels {
      levelPrefixSums.append(UInt32(levelSum))
      levelSum += level.count
    }
    for levelID in levels.indices {
      let prefix = levelPrefixSums[levelID]
      for i in mappedIndices[levelID].indices {
        mappedIndices[levelID][i] += prefix
      }
    }
    
    // Merge the hierarchy levels into a consecutive array.
    nodes = []
    for levelID in levels.indices {
      for var node in levels[levelID] {
        if levelID > 0 {
          node.parentIndex += levelPrefixSums[levelID - 1]
        }
        if !isMax(node.branchesMask) {
          node.branchesIndex += levelPrefixSums[levelID + 1]
        }
        nodes.append(node)
      }
    }
    
    return mappedIndices.flatMap { $0 }
  }
}
