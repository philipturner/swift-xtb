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

/// An octree data structure designed for efficient traversal.
public struct Octree {
  /// The cells from every hierarchy level, in Morton order.
  ///
  /// Each array element is a tuple:
  /// - nextElement: If this element terminates the current 2x2x2 cell, the next
  ///                element is 'nil'. Otherwise, it points to the location
  ///                right after the children and sub-children.
  /// - childCount: Each element is followed by 8 elements specifying its
  ///               children. The number of children must be either 0 or 8.
  public var linkedList: [(nextElement: UInt32?, childCount: UInt8)] = []
  
  /// Center (first three lanes) and grid spacing (fourth lane) of each cell.
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
  
  // Efficient function to add children to several nodes at once.
  mutating func insertChildren(at nodes: [UInt32]) {
    var insertionMarks = [Bool](repeating: false, count: linkedList.count)
    
    // Assert that each specified node exists, and its child count is zero.
    for nodeID in nodes {
      guard nodeID < linkedList.count else {
        fatalError("Element does not exist.")
      }
      let element = linkedList[Int(nodeID)]
      guard element.childCount == 0 else {
        fatalError("Attempted to insert children that already exist.")
      }
      insertionMarks[Int(nodeID)] = true
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
      
      let mark = insertionMarks[currentNodeID]
      if mark == true {
        nextElementOffset += 8
      }
    }
    
    // Create a new array by scanning the current one serially.
    var newLinkedList: [(nextElement: UInt32?, childCount: UInt8)] = []
    var newMetadata: [SIMD4<Float>] = []
    for currentNodeID in linkedList.indices {
      var element = linkedList[currentNodeID]
      let mark = insertionMarks[currentNodeID]
      if mark == true {
        element.childCount = 8
      }
      newLinkedList.append(element)
      if mark == true {
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
      if mark == true {
        for childID in 0..<8 {
          let xIndex = UInt32(childID) % 2
          let yIndex = UInt32(childID >> 1) % 2
          let zIndex = UInt32(childID >> 2) % 2
          var delta: SIMD3<Float> = .init(repeating: -0.25)
          let indices = SIMD3<UInt32>(xIndex, yIndex, zIndex)
          delta.replace(with: 0.25, where: indices .> 0)
          
          let parentPosition = unsafeBitCast(
            parentMetadata, to: SIMD3<Float>.self)
          let parentSpacing = parentMetadata.w
          let childPosition = parentPosition + delta * parentSpacing
          let childSpacing = parentSpacing / 2
          newMetadata.append(SIMD4(childPosition, childSpacing))
        }
      }
    }
    
    // Replace the current arrays with the expanded ones.
    linkedList = newLinkedList
    metadata = newMetadata
  }
  
  // Efficient function to remove children from several nodes at once.
  mutating func removeChildren(from nodes: [UInt32]) {
    var removalMarks = [Bool](repeating: false, count: linkedList.count)
    
    // Assert that each specified node exists, and its child count is 8.
    for nodeID in nodes {
      guard nodeID < linkedList.count else {
        fatalError("Element does not exist.")
      }
      let element = linkedList[Int(nodeID)]
      guard element.childCount == 0 else {
        fatalError("Child cells cannot be removed directly.")
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
      let mark = removalMarks[currentNodeID]
      let pointingNodeID = pointingNodeMap[currentNodeID]
      
      if mark == true {
        nextElementOffset -= 1
      } else if let pointingNodeID {
        let pointingElement = linkedList[Int(pointingNodeID)]
        guard let nextElement = pointingElement.nextElement else {
          fatalError("This should never happen.")
        }
        
        let newValue = nextElement + UInt32(nextElementOffset)
        linkedList[Int(pointingNodeID)].nextElement = newValue
      }
    }
    
    // Remove the array elements for child cells.
    var newLinkedList: [(nextElement: UInt32?, childCount: UInt8)] = []
    var newMetadata: [SIMD4<Float>] = []
    for currentNodeID in linkedList.indices {
      let mark = removalMarks[currentNodeID]
      guard mark == false else {
        continue
      }
      newLinkedList.append(linkedList[currentNodeID])
      newMetadata.append(metadata[currentNodeID])
    }
    
    // Replace the current arrays with the expanded ones.
    linkedList = newLinkedList
    metadata = newMetadata
  }
}
