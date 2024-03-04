//
//  MTLCompiler.swift
//
//
//  Created by Philip Turner on 3/3/24.
//

import Metal
import QuartzCore

public class MTLCompiler {
  public internal(set) var buildProductsDirectory: URL
  var xcodePath: String = "/Applications/Xcode 14.2.app"
  
  public init() {
    let temporaryDirectory = FileManager.default.temporaryDirectory
    buildProductsDirectory = temporaryDirectory
      .appendingPathComponent("metal-compiler")
  }
  
  public func setXcodePath(_ path: String) {
    self.xcodePath = path
  }
  
  // This has about 200 ms latency, so invoke it sparingly. In the final
  // version of this library, we'll compile all GPU shaders beforehand and cache
  // the Metal resource objects at runtime.
  //
  // Near-term, create a key-value cache of the shader sources to the
  // compiled objects. This could speed up some workflows during development of
  // the library.
  public func compile(_ source: String) -> MTLLibrary {
    //
    //  build.swift
    //  MetalFlashAttention
    //
    //  Created by Philip Turner on 6/26/23.
    //
    
    // MARK: - Parse Arguments
    
    struct BuildSettings {
      var externalMetallibPath: String? = nil
      var platform: Platform? = nil
      var verbose: Bool = false
      var xcodePath: String = "/Applications/Xcode 14.2.app"
      
      enum Platform {
        case iOS
        case macOS
        
        var metalToolsPath: String {
          switch self {
          case .iOS:
            return "ios"
          case .macOS:
            return "macos"
          }
        }
        
        var deploymentVersionArgument: String {
          switch self {
          case .iOS:
            return "-mios-version-min=16.0.0"
          case .macOS:
            return "-mmacosx-version-min=13.0.0"
          }
        }
        
        var xcrunSDK: String {
          switch self {
          case .iOS:
            return "iphoneos"
          case .macOS:
            return "macosx"
          }
        }
      }
      
      func metalToolPath(executable: String) -> String {
        guard let metalToolsPath = platform?.metalToolsPath else {
          fatalError("Must specify platform before locating Metal tools.")
        }
        var output = xcodePath
        output += "/Contents/Developer/Toolchains/XcodeDefault.xctoolchain"
        output += "/usr/metal/\(metalToolsPath)/bin/\(executable)"
        return output
      }
      
      func xcrunMetalArguments(executable: String) -> [String] {
        guard let xcrunSDK = platform?.xcrunSDK else {
          fatalError("Must specify platform before locating Metal tools.")
        }
        return ["-sdk", xcrunSDK, executable]
      }
    }
    
    var settings = BuildSettings()
    settings.verbose = Bool.random() ? true : true
    settings.platform = .macOS
    settings.xcodePath = self.xcodePath
    
    // MARK: - Prepare File Directories
    
    func directoryExists(url: URL) -> Bool {
      var isDirectory: ObjCBool = false
      let succeeded = FileManager.default.fileExists(
        atPath: url.path, isDirectory: &isDirectory)
      return succeeded && isDirectory.boolValue
    }
    
    func fileExists(url: URL) -> Bool {
      var isDirectory: ObjCBool = false
      let succeeded = FileManager.default.fileExists(
        atPath: url.path, isDirectory: &isDirectory)
      return succeeded && !isDirectory.boolValue
    }
    
    func assertDirectoryExists(url: URL, line: UInt = #line) {
      guard directoryExists(url: url) else {
        fatalError("""
          Line \(line):
          Directory not found at '\(url.path)'.
          """)
      }
    }
    
    func assertFileExists(url: URL, line: UInt = #line) {
      guard fileExists(url: url) else {
        fatalError("""
          Line \(line):
          File not found at '\(url.path)'.
          """)
      }
    }
    
    func touchDirectory(url: URL) {
      if !directoryExists(url: url) {
        try! FileManager.default.createDirectory(
          at: url, withIntermediateDirectories: false)
      }
      assertDirectoryExists(url: url)
    }
    
    let workDir = self.buildProductsDirectory
    touchDirectory(url: workDir)
    
    let buildDir = workDir.appending(component: "build")
    let libDir = buildDir.appending(component: "lib")
    let srcDir = buildDir.appending(component: "src")
    touchDirectory(url: buildDir)
    touchDirectory(url: libDir)
    touchDirectory(url: srcDir)
    
    func createFile(source: String, name: String) {
      guard let data = source.data(using: .utf8) else {
        fatalError("Could not encode source string as UTF-8.")
      }
      let destinationURL = srcDir.appendingPathComponent(name)
      try! data.write(to: destinationURL)
    }
    createFile(
      source: metalSimdgroupEvent,
      name: "metal_simdgroup_event")
    createFile(
      source: metalSimdgroupMatrixStorage,
      name: "metal_simdgroup_matrix_storage")
    createFile(
      source: source,
      name: "File.metal")
    
    let metalURL = srcDir.appendingPathComponent("File.metal")
    let airURL = buildDir.appendingPathComponent("File.air")
    let airPath = airURL.relativePath
    try? FileManager.default.removeItem(atPath: airPath)
    guard FileManager.default.createFile(
      atPath: airPath, contents: nil) else {
      fatalError("Could not create destination path '\(airPath)'.")
    }
    
    // MARK: - Compile AIR File
    
    // Arguments to invoke the Metal compiler with.
    var arguments: [String] = []
    arguments.append(settings.platform!.deploymentVersionArgument)
    arguments.append("-c")
    
    // Suppress compiler warnings unless the user enters '--verbose'.
    if settings.verbose {
      arguments.append("-Wno-unused-function")
      arguments.append("-Wno-unused-variable")
    } else {
      arguments.append("-w")
    }
    
    let process = Process()
    let toolPath = settings.metalToolPath(executable: "metal")
    process.executableURL = URL(filePath: toolPath)
    process.arguments = arguments + [metalURL.path, "-o", airPath]
    try! process.run()
    process.waitUntilExit()
    if process.terminationStatus != 0 {
      fatalError("Could not compile source.")
    }
    
    // MARK: - Compile Metal Library
    
    arguments = []
    
    // Package the metallib using the up-to-date Xcode version.
    func runProcess() {
      let process = try! Process.run(
        URL(fileURLWithPath: "/usr/bin/xcrun"),
        arguments: arguments)
      process.waitUntilExit()
    }
    
    let metallibName = "File.metallib"
    let metallibURL = libDir.appending(component: metallibName)
    arguments = []
    arguments += settings.xcrunMetalArguments(executable: "metal")
    arguments.append(settings.platform!.deploymentVersionArgument)
    arguments.append(airPath)
    arguments.append("-o")
    arguments.append(metallibURL.path)
    runProcess()
    
    // MARK: - Instantiate Resource Object
    
    let device = MTLCreateSystemDefaultDevice()!
    let library = try! device.makeLibrary(URL: metallibURL)
    
    try! FileManager.default.removeItem(at: workDir)
    return library
  }
}
