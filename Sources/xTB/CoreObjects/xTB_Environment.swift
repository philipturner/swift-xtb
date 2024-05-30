//
//  xTB_Environment.swift
//
//
//  Created by Philip Turner on 5/29/24.
//

/// Calculation environment.
public class xTB_Environment {
  var env: xtb_TEnvironment?
  
  /// Create new xtb calculation environment object.
  public init() {
    env = xtb_newEnvironment()
    guard env != nil else {
      fatalError("Could not create new xTB_Environment.")
    }
  }
  
  /// Delete a xtb calculation environment object.
  deinit {
    xtb_delEnvironment(&env)
  }
  
  /// Check current status of calculation environment.
  public var status: Int {
    let status = xtb_checkEnvironment(env)
    return Int(status)
  }
  
  /// Show and empty error stack.
  public func show() {
    xtb_showEnvironment(env, nil)
  }
  
  /// Possible print levels for API calls.
  public enum Verbosity: UInt32 {
    case full = 2
    case minimal = 1
    case muted = 0
  }
  
  /// Set verbosity of calculation output.
  ///
  /// WARNING: Never call the getter. There is no corresponding xTB function.
  public var verbosity: Verbosity {
    get {
      fatalError("The property 'verbosity' has no getter.")
    }
    set {
      xtb_setVerbosity(env, Int32(newValue.rawValue))
    }
  }
}
