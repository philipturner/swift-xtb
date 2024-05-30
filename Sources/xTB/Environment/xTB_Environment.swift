//
//  xTB_Environment.swift
//
//
//  Created by Philip Turner on 5/29/24.
//

/// Calculation environment.
public class xTB_Environment {
  static let _environment: xtb_TEnvironment = {
    guard let env = xtb_newEnvironment() else {
      fatalError("Could not create new xTB_Environment.")
    }
    return env
  }()
}

extension xTB_Environment {
  /// Possible print levels for API calls.
  public enum Verbosity: UInt32 {
    case full = 2
    case minimal = 1
    case muted = 0
  }
  
  /// Verbosity of calculation output.
  ///
  /// The default value is `.full`.
  public static var verbosity: Verbosity = .full {
    didSet {
      print("Changing verbosity to \(verbosity)")
      xtb_setVerbosity(_environment, Int32(verbosity.rawValue))
    }
  }
  
  /// Check current status of calculation environment.
  public static var status: Int {
    let status = xtb_checkEnvironment(_environment)
    return Int(status)
  }
  
  /// Show and empty error stack.
  public static func show() {
    xtb_showEnvironment(_environment, nil)
  }
}
