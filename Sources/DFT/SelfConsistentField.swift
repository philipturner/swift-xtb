//
//  SelfConsistentField.swift
//
//
//  Created by Philip Turner on 2/4/24.
//

import OpenCL

public class SelfConsistentField {
  @usableFromInline
  var _queue: CLCommandQueue?
  
  @usableFromInline
  var _ram: RAM?
  
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
