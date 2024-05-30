//
//  xTB_Environment.swift
//
//
//  Created by Philip Turner on 5/29/24.
//

/// Calculation environment.
public class xTB_Environment {
  var pointer: xtb_TEnvironment
  
  /// Possible print levels for API calls.
  public enum Verbosity: UInt32 {
    case full = 2
    case minimal = 1
    case muted = 0
  }
  
  var _verbosity: Verbosity = .full
  
  /// Create new xtb calculation environment object.
  public init() {
    guard let env = xtb_newEnvironment() else {
      fatalError("Could not create new xTB_Environment.")
    }
    pointer = env
  }
  
  /// Delete a xtb calculation environment object.
  deinit {
    xtb_delEnvironment(&pointer)
  }
}

extension xTB_Environment {
  /// Check current status of calculation environment.
  public var status: Int {
    let status = xtb_checkEnvironment(pointer)
    return Int(status)
  }
  
  /// Show and empty error stack.
  public func show() {
    xtb_showEnvironment(pointer, nil)
  }
  
  /// Verbosity of calculation output.
  ///
  /// The default value is `.full`.
  public var verbosity: Verbosity {
    get {
      _verbosity
    }
    set {
      _verbosity = newValue
      xtb_setVerbosity(pointer, Int32(newValue.rawValue))
    }
  }
}
