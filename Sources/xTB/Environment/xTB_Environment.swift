//
//  xTB_Environment.swift
//  swift-xtb
//
//  Created by Philip Turner on 5/29/24.
//

import C_xTB

/// Calculation environment.
public class xTB_Environment {
  // Lazily initialized global variable.
  static let tEnvironment: xtb_TEnvironment = {
    return xTB_Environment.createObject()
  }()
  
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
    // TODO: Expose xtb_getError for cleaner API to handle errors.
    
    
//    if status != 0 {
//      var characters: [CChar] = []
//      var iterationCount: Int = 0
//      while true {
//        iterationCount += 1
//        print("getError iteration count: \(iterationCount)")
//        
//        var buffer = [CChar](repeating: 0, count: 512)
//        var bufferSize: Int32 = 512
//        xtb_getError(tEnvironment, &buffer, &bufferSize)
//        print("bufferSize:", bufferSize)
//        
//        characters += buffer
//        
//        if buffer[511] == 0 {
//          break
//        }
//      }
//      print("total character count:", characters.count)
//      print("message:", String(cString: characters))
      
//    }
    
    xtb_showEnvironment(tEnvironment, nil)
  }
  
  // Redirect the output to something other than the console.
  //
  // `setOutput` keeps appending to a list of output files. At the end of the program, all of them get written to. Unless the list has 1 element, and that element is `/dev/null`. In that case, the console output is permanently disabled for GFN2-xTB. And conditionally disabled for GFN-FF (if the verbosity also happens to be `.muted`).
  //
  // `releaseOutput` doesn't have any effects on file writing behavior. It doesn't make the console available again.
  public static func setOutput(_ filename: String) {
    xtb_setOutput(tEnvironment, filename)
  }
}
