//
//  ReflectorGeneration.swift
//
//
//  Created by Philip Turner on 4/6/24.
//

import Foundation

// A configuration for a reflector generation.
struct ReflectorGenerationDescriptor {
  var source: UnsafePointer<Float>?
  var dimension: Int?
}

// Creates a reflector that transforms the provided vector into [1, 0, 0, ...].
struct ReflectorGeneration {
  // TODO: Try to write the reflector directly into a user-specified memory
  // allocation. Attempt this after refactoring 'BulgeChasing' to use the new
  // API.
  var reflector: [Float]
  
  // TODO: Try to fuse tau directly into vector, instead of leaving it as a
  // separate variable. Attempt this after building the 'WYTransform' API.
  var tau: Float = .zero
  
  init(descriptor: ReflectorGenerationDescriptor) {
    guard let source = descriptor.source,
          let dimension = descriptor.dimension else {
      fatalError("Descriptor was incomplete.")
    }
    reflector = Array(repeating: .zero, count: dimension)
    
    // Take the norm of the vector.
    var norm: Float = .zero
    for elementID in 0..<dimension {
      norm += source[elementID] * source[elementID]
    }
    norm.formSquareRoot()
    
    // Check for NANs.
    let oldSubdiagonal = source[0]
    let newSubdiagonal = norm * Float((oldSubdiagonal >= 0) ? -1 : 1)
    let epsilon: Float = 2 * .leastNormalMagnitude
    guard (newSubdiagonal - oldSubdiagonal).magnitude > epsilon else {
      return
    }
    guard newSubdiagonal.magnitude > epsilon else {
      return
    }
    
    // Predict the normalization factor.
    tau = (newSubdiagonal - oldSubdiagonal) / newSubdiagonal
    let normalizationFactor = 1 / (oldSubdiagonal - newSubdiagonal)
    
    // Write to the reflector.
    for elementID in 0..<dimension {
      let element = source[elementID]
      reflector[elementID] = element * normalizationFactor
    }
    reflector[0] = 1
  }
}
