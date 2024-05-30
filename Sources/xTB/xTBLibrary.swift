//
//  xTBLibrary.swift
//  
//
//  Created by Philip Turner on 5/29/24.
//

#if canImport(Darwin)
import Darwin
#elseif canImport(Glibc)
import Glibc
#elseif os(Windows)
import CRT
import WinSDK
#endif

//===----------------------------------------------------------------------===//
// The `xTBLibrary` struct that loads xTB symbols at runtime.
//===----------------------------------------------------------------------===//

public struct xTBLibrary {
  public enum Error: Swift.Error, Equatable, CustomStringConvertible {
    case xtbLibraryNotFound
    
    public var description: String {
      switch self {
      case .xtbLibraryNotFound:
        return "xTB library not found. Set the library path to an xTB library."
      }
    }
  }
  
  private static var isXTBLibraryLoaded = false
  private static var _xtbLibraryHandle: UnsafeMutableRawPointer?
  private static var xtbLibraryHandle: UnsafeMutableRawPointer? {
    try! xTBLibrary.loadLibrary()
    return self._xtbLibraryHandle
  }
  
  /// Tries to load the xTB library, will throw an error if no compatible
  /// library is found.
  public static func loadLibrary() throws {
    guard !isXTBLibraryLoaded else { return }
    let xtbLibraryHandle = self.loadXTBLibrary()
    guard self.isXTBLibraryLoaded(at: xtbLibraryHandle) else {
      throw Error.xtbLibraryNotFound
    }
    self.isXTBLibraryLoaded = true
    self._xtbLibraryHandle = xtbLibraryHandle
  }
  
  public static func loadSymbol<T>(
    name: String,
    type: T.Type = T.self
  ) -> T {
    print("Loading symbol '\(name)' from the xTB library...")
    return unsafeBitCast(
      self.loadSymbol(self.xtbLibraryHandle, name), to: type)
  }
}

// Methods of `xTBLibrary` required to load the xTB library.
extension xTBLibrary {
  private static var libraryPath: String?
  
  private static func loadSymbol(
    _ libraryHandle: UnsafeMutableRawPointer?,
    _ name: String
  ) -> UnsafeMutableRawPointer? {
#if canImport(Darwin) || canImport(Glibc)
    return dlsym(libraryHandle, name)
#elseif os(Windows)
    guard let libraryHandle = libraryHandle else { return nil }
    let moduleHandle = libraryHandle
      .assumingMemoryBound(to: HINSTANCE__.self)
    let moduleSymbol = GetProcAddress(moduleHandle, name)
    return unsafeBitCast(moduleSymbol, to: UnsafeMutableRawPointer?.self)
#endif
  }
  
  private static func isXTBLibraryLoaded(
    at xtbLibraryHandle: UnsafeMutableRawPointer? = nil
  ) -> Bool {
    return self.loadSymbol(xtbLibraryHandle, "xtb_getAPIVersion") != nil
  }
  
  private static func loadXTBLibrary() -> UnsafeMutableRawPointer? {
    var xtbLibraryHandle: UnsafeMutableRawPointer?
    if let xtbLibraryPath = xTBLibrary.libraryPath {
      xtbLibraryHandle = self.loadXTBLibrary(at: xtbLibraryPath)
    }
    return xtbLibraryHandle
  }
  
  private static func loadXTBLibrary(
    at path: String
  ) -> UnsafeMutableRawPointer? {
    print("Trying to load library at '\(path)'...")
#if canImport(Darwin) || canImport(Glibc)
    let xtbLibraryHandle = dlopen(path, RTLD_LAZY | RTLD_GLOBAL)
#elseif os(Windows)
    let xtbLibraryHandle = UnsafeMutableRawPointer(LoadLibraryA(path))
#endif
    if xtbLibraryHandle != nil {
      print("Library at '\(path)' was successfully loaded.")
    } else {
      print("Library at '\(path)' failed to load.")
    }
    return xtbLibraryHandle
  }
}

// Methods of `xTBLibrary` required to set a given xTB version or library path.
extension xTBLibrary {
  private static func enforceNonLoadedXTBLibrary(
    function: String = #function
  ) {
    guard !self.isXTBLibraryLoaded else {
      fatalError("""
        Error: \(function) should not be called after any xTB library \
        has already been loaded.
        """)
    }
  }
  
  // Use the xTB library at the specified path.
  // - Parameter path: Path of the xTB library to load or nil to use the
  //   default search path.
  public static func useLibrary(at path: String?) {
    self.enforceNonLoadedXTBLibrary()
    xTBLibrary.libraryPath = path
  }
}
