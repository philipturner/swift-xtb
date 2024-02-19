//
//  Octree.swift
//
//
//  Created by Philip Turner on 2/18/24.
//

// An octree data structure optimized for traversal on highly parallel computer
// architectures (CPU and GPU).
//
// Stores a set of fragments. Each fragment internally stores a 2x2x2 group
// of sub-fragments contiguously. One can compute the first erivative at the
// center of the cell for XC functionals.

struct Octree {
  // Stores the child cells' data. Even if there aren't child nodes, there are
  // always child grid cells.
  var cellValues: [SIMD8<Float>] = []
  
  // If this element terminates the current 2x2x2 cell, the next element
  // is 'nil'. Otherwise, it points to the location right after the children
  // and sub-children.
  //
  // Each element is followed by 8 elements specifying its children. The number
  // of children must be either 0 or 8.
  var linkedList: [(nextElement: UInt32?, childCount: UInt8)] = []
  
  init() {
    cellValues = [SIMD8.zero]
    linkedList = [(nextElement: nil, childCount: 0)]
  }
  
  // Efficient function to add children to several nodes at once.
  mutating func appendChildren(to nodes: [UInt32]) {
    var insertionMarks = [Bool](repeating: false, count: linkedList.count)
    
    // Assert that each specified node exists, and its child count is zero.
    for nodeID in nodes {
      guard nodeID < linkedList.count else {
        fatalError("Element does not exist.")
      }
      let element = linkedList[Int(nodeID)]
      guard element.childCount == 0 else {
        // Converse: Removed children from vacant cell.
        fatalError("Inserted children into occupied cell.")
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
    var newCellValues: [SIMD8<Float>] = []
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
      
      // Initialize the child cells' data to SIMD8<Float>.zero. Do not modify
      // the data currently belonging to the expanded cell.
      let value = cellValues[Int(currentNodeID)]
      newCellValues.append(value)
      if mark == true {
        for _ in 0..<8 {
          let value: SIMD8<Float> = .zero
          newCellValues.append(value)
        }
      }
    }
    
    // Replace the current arrays with the expanded ones.
    cellValues = newCellValues
    linkedList = newLinkedList
  }
  
  // Efficient function to remove children from several nodes at once.
  mutating func removeChildren(to nodes: [UInt32]) {
    var removalMarks = [Bool](repeating: false, count: linkedList.count)
    
    // Assert that each specified node exists, and its child count is 8.
    for nodeID in nodes {
      guard nodeID < linkedList.count else {
        fatalError("Element does not exist.")
      }
      let element = linkedList[Int(nodeID)]
      guard element.childCount == 0 else {
        fatalError("Removed children from vacant cell.")
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
    
    // Remove the array elements for child cells' data.
    var newLinkedList: [(nextElement: UInt32?, childCount: UInt8)] = []
    var newCellValues: [SIMD8<Float>] = []
    for currentNodeID in linkedList.indices {
      let mark = removalMarks[currentNodeID]
      guard mark == false else {
        continue
      }
      let element = linkedList[currentNodeID]
      let value = cellValues[currentNodeID]
      newLinkedList.append(element)
      newCellValues.append(value)
    }
    
    // Replace the current arrays with the expanded ones.
    cellValues = newCellValues
    linkedList = newLinkedList
  }
}
