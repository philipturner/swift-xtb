//
//  Octree.swift
//
//
//  Created by Philip Turner on 2/18/24.
//

/// A configuration for an octree.
struct OctreeDescriptor {
  /// Required. The power-2 size of the coarsest level.
  var sizeExponent: Int?
}

// The data needed for locating this node in 3D space and memory.
public struct OctreeNode {
  public var centerAndSpacing: SIMD4<Float>
  
  // The average 3D position of the cell's contents.
  @_transparent public var center: SIMD3<Float> {
    unsafeBitCast(centerAndSpacing, to: SIMD3<Float>.self)
  }
  
  // The cube root of the volume.
  @_transparent public var spacing: Float {
    centerAndSpacing.w
  }
  
  // The number of additional elements between this node and its neighbors.
  public var childrenBefore: UInt32
  public var childrenAfter: UInt32
  
  // The index in the parent node. Consecutive nodes in memory may not have
  // contiguous parent indices, because some nodes aren't expanded. Use this to
  // compute the parent's center, then find the first node with that center.
  public var indexInParent: UInt8
  
  // If any children have grand-children, they are marked in this mask. The
  // next array elements will be the child nodes in compacted order. If the
  // child is a leaf node, it isn't marked.
  public var branchesMask: UInt8
}

/// An octree data structure designed for efficient traversal.
public struct Octree {
  /// The cells from every hierarchy level, in Morton order.
  public internal(set) var nodes: [OctreeNode] = []
  
  init(descriptor: OctreeDescriptor) {
    guard let sizeExponent = descriptor.sizeExponent else {
      fatalError("Descriptor was invalid.")
    }
    
    let origin: SIMD3<Float> = .zero
    let size = Float(sign: .plus, exponent: sizeExponent, significand: 1)
    let node = OctreeNode(
      centerAndSpacing: SIMD4(origin, size),
      childrenBefore: 0,
      childrenAfter: 0,
      indexInParent: 0,
      branchesMask: 0b0000_0000)
    nodes = [node]
  }
  
  // Efficient function to expand/contract several nodes at once.
  // - expand: The nodes that should acquire children. Each array element is
  //           the parent address and the mask of children that should be
  //           expanded.
  // - contract: The parents who should no longer have children. The children
  //             cannot have grandchildren.
  // - returns: The old nodes' positions in the new list. If no such position
  //            exists, the array element is `UInt32.max`.
  mutating func resizeNodes(
    expand: [(UInt32, UInt8)], contract: [UInt32]
  ) -> [UInt32] {
    var insertionMarks = [UInt8](repeating: 0b0000_0000, count: nodes.count)
    var removalMarks = [Bool](repeating: false, count: nodes.count)
    
    // Assert that each specified node exists, and its child count is zero.
    for (nodeID, expansionMask) in expand {
      guard nodes.indices.contains(Int(nodeID)) else {
        fatalError("Parent does not exist.")
      }
      let node = nodes[Int(nodeID)]
      guard expansionMask & node.branchesMask == 0 else {
        fatalError("Attempted to insert children that already exist.")
      }
      insertionMarks[Int(nodeID)] = expansionMask
    }
    
    // Assert that each specified node exists, and its children don't have
    // grandchildren.
    for nodeID in contract {
      guard nodes.indices.contains(Int(nodeID)) else {
        fatalError("Parent does not exist.")
      }
      let node = nodes[Int(nodeID)]
      guard node.branchesMask == 0 else {
        fatalError("Attempted to contract node with grandchildren.")
      }
      guard insertionMarks[Int(nodeID)] == 0 else {
        fatalError("Attempted to expand and contract a node at the same time.")
      }
      removalMarks[Int(nodeID)] = true
    }
    
    // Modify the existing linked list so 'nextElement' already points to
    // where it should in the new list.
    //
    // TODO: Something similar, but tracking the number of elements expanded or
    // contracted between a specific node pair. This could get very complicated
    // if we want it to be an in-place operation.
    //
    // Alternatively, try updating from left to right. Store a stack of the last
    // time a particular level of the hierarchy was hit. That instantly
    // results in a O(n) childrenBefore computation. Repeat a similar process
    // for the other traversal direction for childrenAfter.
    // - WARNING: Need to erase certain nodes of the stack once they go out of
    //            context (belong to the wrong parent). I think this happens
    //            automatically when popping the stack to reach a higher level.
    var nextElementOffset: Int = .zero
    for currentNodeID in linkedList.indices {
      let pointingNodeID = pointingNodeMap[currentNodeID]
      if let pointingNodeID {
        let pointingElement = linkedList[Int(pointingNodeID)]
        guard let nextElement = pointingElement.nextElement else {
          fatalError("This should never happen.")
        }
        
        let newValue = nextElement + UInt32(nextElementOffset)
        linkedList[Int(pointingNodeID)].nextElement = newValue
      }
      
      if insertionMarks[currentNodeID] {
        nextElementOffset += 8
      }
      if removalMarks[currentNodeID] {
        nextElementOffset -= 1
      }
    }
    
    // Create a new array by scanning the current one serially.
    var mappedPositions = [UInt32](repeating: .max, count: linkedList.count)
    var newLinkedList: [(nextElement: UInt32?, childCount: UInt8)] = []
    var newMetadata: [SIMD4<Float>] = []
    for currentNodeID in linkedList.indices {
      var element = linkedList[currentNodeID]
      if insertionMarks[currentNodeID] {
        element.childCount = 8
      }
      if removalMarks[currentNodeID] {
        continue
      }
      if !insertionMarks[currentNodeID] {
        mappedPositions[currentNodeID] = UInt32(newLinkedList.count)
      }
      newLinkedList.append(element)
      if insertionMarks[currentNodeID] {
        for childID in 0..<8 {
          var nextElement: Optional = UInt32(newLinkedList.count)
          if childID == 7 {
            nextElement = nil
          }
          newLinkedList.append((nextElement: nextElement, childCount: 0))
        }
      }
      
      let parentMetadata = metadata[Int(currentNodeID)]
      newMetadata.append(parentMetadata)
      if insertionMarks[currentNodeID] {
        var x = SIMD8<Float>(0, 1, 0, 1, 0, 1, 0, 1) * 0.5 - 0.25
        var y = SIMD8<Float>(0, 0, 1, 1, 0, 0, 1, 1) * 0.5 - 0.25
        var z = SIMD8<Float>(0, 0, 0, 0, 1, 1, 1, 1) * 0.5 - 0.25
        x = x * parentMetadata.w + parentMetadata.x
        y = y * parentMetadata.w + parentMetadata.y
        z = z * parentMetadata.w + parentMetadata.z
        
        for childID in 0..<8 {
          let parentSpacing = parentMetadata.w
          let childPosition = SIMD3(x[childID], y[childID], z[childID])
          let childSpacing = parentSpacing / 2
          newMetadata.append(SIMD4(childPosition, childSpacing))
        }
      }
    }
    
    // Replace the current arrays with the expanded ones.
    linkedList = newLinkedList
    metadata = newMetadata
    return mappedPositions
  }
}
