//
//  Octree.swift
//
//
//  Created by Philip Turner on 2/18/24.
//

import Foundation

// An octree data structure optimized for traversal on highly parallel computer
// architectures (CPU and GPU).
//
// Stores a set of fragments. Each fragment internally stores a 2x2x2 group
// of sub-fragments contiguously. It can compute the first derivative at the
// center of the cell for XC functionals.

struct OctreeDescriptor {
  // The minimum bounding volume of the highest octree level. This number can
  // be very liberal.
  var lowerCorner: SIMD3<Float>?
  var upperCorner: SIMD3<Float>?
}

struct Octree {
  init(descriptor: OctreeDescriptor) {
    
  }
}
