//
//  xTB_Environment.swift
//  swift-xtb
//
//  Created by Philip Turner on 5/29/24.
//

/// Calculation environment.
public class xTB_Environment {
  static let tEnvironment: xtb_TEnvironment = {
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
  /// The default value is `.minimal`.
  public static var verbosity: Verbosity = .minimal {
    didSet {
      xtb_setVerbosity(tEnvironment, Int32(verbosity.rawValue))
    }
  }
  
  /// Check current status of calculation environment.
  public static var status: Int {
    let status = xtb_checkEnvironment(tEnvironment)
    return Int(status)
  }
  
  /// Show and empty error stack.
  public static func show() {
    xtb_showEnvironment(tEnvironment, nil)
  }
  
  /// Redirect the output to something other than stdout.
  ///
  /// `setOutput` keeps appending to a list of output files. At the end of the
  /// program, all of them get written to. Unless the file happens to be
  /// `/dev/null`. In that case, the output is permanently disabled for
  /// GFN2-xTB. And conditionally disabled for GFN-FF (if the verbosity also
  /// happens to be `.muted`).
  ///
  /// `releaseOutput` doesn't do anything, at least on macOS.
  public static func setOutput(_ filename: String) {
    xtb_setOutput(tEnvironment, filename)
  }
}
