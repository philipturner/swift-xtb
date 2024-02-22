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

struct OctreeNode {
  // The data in the octree. For vector lanes corresponding to branches, this is
  // implicitly assumed to be the weighted average of all children. If the
  // assumption is not correct, you must account for that in calling code.
  var data: SIMD8<Float>
  var centerAndSpacing: SIMD4<Float>
  
  // The average 3D position of the cell's contents.
  @_transparent var center: SIMD3<Float> {
    unsafeBitCast(centerAndSpacing, to: SIMD3<Float>.self)
  }
  
  // The cube root of the volume.
  @_transparent var spacing: Float {
    centerAndSpacing.w
  }
  
  // A means for traversing to neighboring nodes.
  var childrenBefore: UInt32
  var childrenAfter: UInt32
  
  // The index in the parent node. Consecutive nodes in memory may not have
  // contiguous parent indices, because some nodes aren't expanded.
  var indexInParent: UInt8
  
  // If any children have sub-children, they are marked in this mask. The
  // consecutive array elements contain the child nodes in compacted order.
  //
  // Otherwise, the child is implicitly specified to exist, but it terminates
  // the tree. Places where the mask evaluates to 0 are "sources of truth".
  var branchesMask: UInt8
}

/// An octree data structure designed for efficient traversal.
public struct Octree {
  /// The cells from every hierarchy level, in Morton order.
  ///
  /// Each array element is a tuple:
  /// - nextElement: If this element terminates the current 2x2x2 cell, the next
  ///                element is 'nil'. Otherwise, it points to the location
  ///                right after the children.
  /// - previousElement: If this element begins the current 2x2x2 cell, the
  ///                    previous element is 'nil'. Otherwise, it points to the
  ///                    location right before the previous node's children.
  /// - branches: If any children have sub-children, they are marked in this
  ///              mask. The consecutive array elements contain the child nodes
  ///              in compacted order. Otherwise, the child is implicitly
  ///              specified to exist, but it terminates the tree. Places where
  ///              the mask evaluates to 0 are "sources of truth".
  public var linkedList: [(nextElement: UInt32?, branches: UInt8)] = []
  
  /// Center (first three lanes) and grid spacing (fourth lane) of the parent
  /// cell.
  public var metadata: [SIMD4<Float>] = []
  
  init(descriptor: OctreeDescriptor) {
    guard let sizeExponent = descriptor.sizeExponent else {
      fatalError("Descriptor was invalid.")
    }
    
    let origin: SIMD3<Float> = .zero
    let size = Float(sign: .plus, exponent: sizeExponent, significand: 1)
    metadata = [SIMD4(origin, size)]
    linkedList = [(nextElement: nil, childCount: 0)]
  }
  
  // Efficient function to expand/contract several modes at once.
  //
  // Returns the old nodes' positions in the new list. If no such position
  // exists, the array element is `UInt32.max`.
  mutating func modifyNodes(expand: [UInt32], contract: [UInt32]) -> [UInt32] {
    var insertionMarks = [Bool](repeating: false, count: linkedList.count)
    var removalMarks = [Bool](repeating: false, count: linkedList.count)
    
    // Assert that each specified node exists, and its child count is zero.
    for nodeID in expand {
      guard nodeID < linkedList.count else {
        fatalError("Element does not exist.")
      }
      let element = linkedList[Int(nodeID)]
      guard element.childCount == 0 else {
        fatalError("Attempted to insert children that already exist: \(expand).")
      }
      insertionMarks[Int(nodeID)] = true
    }
    
    // Assert that each specified node exists, and its child count is 8.
    for nodeID in contract {
      guard nodeID < linkedList.count else {
        fatalError("Element does not exist.")
      }
      let element = linkedList[Int(nodeID)]
      guard element.childCount == 8 else {
        fatalError("Contracted cell must have children.")
      }
      linkedList[Int(nodeID)].childCount = 0
      
      // Assert that none of the children have sub-children.
      for childID in 0..<8 {
        let childNodeID = Int(nodeID) + 1 + childID
        let childElement = linkedList[childNodeID]
        guard childElement.childCount == 0 else {
          fatalError("Removed an occupied child cell.")
        }
        removalMarks[childNodeID] = true
      }
    }
    
    // Create a list of the previous elements that reference the current one.
    var pointingNodeMap = [UInt32?](repeating: nil, count: linkedList.count)
    for currentNodeID in linkedList.indices {
      let element = linkedList[currentNodeID]
      if let nextElement = element.nextElement {
        guard pointingNodeMap[Int(nextElement)] == nil else {
          fatalError("This should never happen.")
        }
        pointingNodeMap[Int(nextElement)] = UInt32(currentNodeID)
      }
    }
    
    // Modify the existing linked list so 'nextElement' already points to
    // where it should in the new list.
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
