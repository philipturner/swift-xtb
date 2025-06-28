//
//  xTB_Environment.swift
//  swift-xtb
//
//  Created by Philip Turner on 5/29/24.
//

import C_xTB

public class xTB_Environment {
  /// Lazily initialized singleton for the environment.
  nonisolated(unsafe)
  static let tEnvironment: xtb_TEnvironment = {
    return xTB_Environment.createObject()
  }()
  
  /// Create the reference-counted object from the C API.
  static func createObject() -> xtb_TEnvironment {
    let env = xtb_newEnvironment()
    guard let env else {
      fatalError("Could not create new xTB_Environment.")
    }
    return env
  }
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
  nonisolated(unsafe)
  public static var verbosity: Verbosity = .minimal {
    didSet {
      xtb_setVerbosity(tEnvironment, Int32(verbosity.rawValue))
    }
  }
  
  /// Check for a nonzero status that indicates an error.
  public static var status: Int {
    let status = xtb_checkEnvironment(tEnvironment)
    return Int(status)
  }
  
  /// Show errors to the console and purge the errors.
  public static func show() {
    xtb_showEnvironment(tEnvironment, nil)
  }
  
  /// More reliable method for reading errors, in situations where the console
  /// is disabled by unexpected C API behavior.
  public static func flushErrorStack() -> String {
    // Once you call into the function, all remaining chunks are forfeited.
    // So make sure everything gets retrieved in the first call.
    let BUFFER_CHUNK_SIZE: Int = 4096
    
    var buffer = [UInt8](repeating: 0, count: BUFFER_CHUNK_SIZE)
    var bufferSize = Int32(BUFFER_CHUNK_SIZE)
    xtb_getError(tEnvironment, &buffer, &bufferSize)
    
    guard buffer[BUFFER_CHUNK_SIZE - 1] == 0 else {
      fatalError("Buffer chunk was not null-terminated.")
    }
    return String(decoding: buffer, as: UTF8.self)
  }
  
  /// Redirect the output to something other than the console.
  ///
  /// `setOutput` keeps appending to a list of output files. At the end of the
  /// program, all of them get written to. Unless the list has 1 element, and
  /// that element is `/dev/null`. In that case, the console output is
  /// permanently disabled for GFN2-xTB. And conditionally disabled for GFN-FF
  /// (if the verbosity also happens to be `.muted`).
  ///
  /// `releaseOutput` doesn't have any effects on file writing behavior. It
  /// doesn't make the console available again.
  public static func setOutput(_ filename: String) {
    xtb_setOutput(tEnvironment, filename)
  }
}
