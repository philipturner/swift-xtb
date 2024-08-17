//
//  SelfConsistentField.swift
//
//
//  Created by Philip Turner on 2/4/24.
//

import OpenCL

// TODO: Rewrite all of this. Create an MVP that uses the CPU, then port it to
// the GPU. Ensure the algorithms used are GPU-friendly, even if that harms
// CPU performance.
//
// Archive the current idea and associated tests in the documentation folder.

public class SelfConsistentField {
  @usableFromInline
  var _queue: CLCommandQueue?
  
  @usableFromInline
  var _ram: RAM?
  
  /// The number of statically allocated GPU threads to use for computations.
  /// This should be slightly over the GPU's occupancy.
  public var gpuThreadCount: Int = 0
  
  public init() {
    
  }
}

extension SelfConsistentField {
  public var queue: CLCommandQueue {
    @_transparent
    get {
      return _queue!
    }
    @_transparent
    set {
      _queue = newValue
    }
  }
  
  public var ram: RAM {
    @_transparent
    get {
      return _ram!
    }
    @_transparent
    set {
      _ram = newValue
    }
  }
}
