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
// of sub-fragments contiguously. It can compute the first derivative at the
// center of the cell for XC functionals.

struct OctreeDescriptor {
  // The radius of the highest octree level, in powers of 2. All octrees are
  // centered at the origin.
  var highestLevelSize: Int
}

struct Octree {
  // Stores the child cells' data. Even if there aren't child nodes, there are
  // always child grid cells.
  var cellValues: [SIMD8<Float>] = []
  
  var highestLevelSize: Int
  
  // If this element terminates the current 2x2x2 cell, the next element
  // is 'nil'. Otherwise, it points to the location right after the children
  // and sub-children.
  //
  // Each element is followed by 8 elements specifying its children. The number
  // of children must be either 0 or 8.
  var linkedList: [(nextElement: UInt32?, childCount: UInt8)] = []
  
  init(descriptor: OctreeDescriptor) {
    cellValues = [SIMD8.zero]
    highestLevelSize = descriptor.highestLevelSize
    linkedList = [(nextElement: nil, childCount: 0)]
  }
  
  // Efficient function to add children to several nodes at once.
  mutating func appendChildren(to nodes: [UInt32]) {
    // Assert that each specified node exists, and its child count is zero.
    
    // Initialize the child cells' data to SIMD8<Float>.zero.
  }
  
  // Efficient function to remove children from several nodes at once.
  mutating func removeChildren(to nodes: [UInt32]) {
    // Assert that each specified node exists, and its child count is 8.
    
    // Assert that none of the children have sub-children.
    
    // Remove the array elements for child cells' data.
  }
}
