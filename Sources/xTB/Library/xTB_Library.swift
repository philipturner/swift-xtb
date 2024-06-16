//
//  xTB_Library.swift
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
// The `xTB_Library` struct that loads xTB symbols at runtime.
//===----------------------------------------------------------------------===//

public struct xTB_Library {
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
    try! xTB_Library.loadLibrary()
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
    return unsafeBitCast(
      self.loadSymbol(self.xtbLibraryHandle, name), to: type)
  }
}

// Methods of `xTB_Library` required to load the xTB library.
extension xTB_Library {
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
    if let xtbLibraryPath = xTB_Library.libraryPath {
      xtbLibraryHandle = self.loadXTBLibrary(at: xtbLibraryPath)
    }
    return xtbLibraryHandle
  }
  
  private static func loadXTBLibrary(
    at path: String
  ) -> UnsafeMutableRawPointer? {
#if canImport(Darwin) || canImport(Glibc)
    let xtbLibraryHandle = dlopen(path, RTLD_LAZY | RTLD_GLOBAL)
#elseif os(Windows)
    let xtbLibraryHandle = UnsafeMutableRawPointer(LoadLibraryA(path))
#endif
    return xtbLibraryHandle
  }
}

// Methods of `xTB_Library` required to set a given xTB version or library path.
extension xTB_Library {
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
    xTB_Library.libraryPath = path
  }
}
