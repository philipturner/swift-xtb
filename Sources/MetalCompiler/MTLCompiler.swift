//
//  MTLCompiler.swift
//
//
//  Created by Philip Turner on 3/3/24.
//

import Metal

public class MTLCompiler {
  public internal(set) var buildProductsDirectory: URL
  var xcodePath: String = "/Applications/Xcode 14.2.app"
  
  public init() {
    let temporaryDirectory = FileManager.default.temporaryDirectory
    buildProductsDirectory = temporaryDirectory
      .appendingPathExtension("metal-build")
  }
  
  public func setXcodePath(_ path: String) {
    self.xcodePath = path
  }
  
  public func compile(_ source: String) -> MTLLibrary {
    // For now, we just delegate to the Metal JIT compiler. Once we establish a
    // vector addition test, we can move on to supporting async copies.
    let device = MTLCreateSystemDefaultDevice()!
    let options = MTLCompileOptions()
    options.fastMathEnabled = true // We aren't doing FP64 emulation yet.
    let library = try! device.makeLibrary(source: source, options: options)
    return library
  }
}
