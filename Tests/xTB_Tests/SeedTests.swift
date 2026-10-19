import Foundation
#if os(Windows)
import WinSDK
#endif
import XCTest
import xTB

final class SeedTests: XCTestCase {
  static func createCalculator() -> xTB_Calculator {
    let system: [SIMD4<Float>] = diamondSystem122
    var calculatorDesc = xTB_CalculatorDescriptor()
    calculatorDesc.atomicNumbers = system.map { UInt8($0.w) }
    calculatorDesc.positions = system.map {
      SIMD3($0.x, $0.y, $0.z)
    }
    calculatorDesc.hamiltonian = .tightBinding
    return xTB_Calculator(descriptor: calculatorDesc)
  }

  /// Counts the SCC iterations reported by xTB's console output, by
  /// redirecting stdout to a file while `body` runs. `body` must compute
  /// exactly one singlepoint.
  static func countIterations(_ body: () throws -> Void) throws -> Int {
    let path = NSTemporaryDirectory() + "SeedTests-\(UUID()).txt"
    _ = FileManager.default.createFile(atPath: path, contents: nil)
    defer { try? FileManager.default.removeItem(atPath: path) }
    
    let redirection = StdoutRedirection(path: path)
    xTB_Environment.verbosity = .full
    defer {
      xTB_Environment.verbosity = .muted
      redirection.restore()
    }
    try body()
    
    // "*** convergence criteria satisfied after 13 iterations ***"
    let text = try String(contentsOfFile: path, encoding: .utf8)
    let iterations: [Int] = text.split(separator: "\n").compactMap { line in
      guard line.contains("convergence criteria satisfied after") else {
        return nil
      }
      let words = line.split(separator: " ")
      guard let index = words.firstIndex(of: "after") else {
        return nil
      }
      return Int(words[index + 1])
    }
    XCTAssertEqual(iterations.count, 1, "Expected exactly one singlepoint.")
    return try XCTUnwrap(iterations.first)
  }
  
  func testSeedAvailability() throws {
    xTB_Environment.verbosity = .muted
    let calculator = SeedTests.createCalculator()
    
    // No singlepoint has been computed yet.
    XCTAssertNil(calculator.seed)
    
    // The singlepoint is up to date.
    _ = calculator.energy
    XCTAssertNotNil(calculator.seed)
    
    // The singlepoint became stale.
    calculator.molecule.positions[0].x += 0.001
    XCTAssertNil(calculator.seed)
  }
  
  func testWarmStartEnergy() throws {
    xTB_Environment.verbosity = .muted
    let warmCalculator = SeedTests.createCalculator()
    let coldCalculator = SeedTests.createCalculator()
    
    _ = warmCalculator.energy
    let seed = try XCTUnwrap(warmCalculator.seed)
    
    // Warm start several displacements from the same seed, like a Hessian.
    let originalPositions = warmCalculator.molecule.positions
    let displacements: [Float] = [0.001, -0.001]
    let warmEnergies = try warmCalculator.withSeed(seed) {
      displacements.map { displacement in
        warmCalculator.molecule.positions = originalPositions
        warmCalculator.molecule.positions[0].x += displacement
        return warmCalculator.energy
      }
    }
    
    for (displacement, warmEnergy) in zip(displacements, warmEnergies) {
      coldCalculator.molecule.positions = originalPositions
      coldCalculator.molecule.positions[0].x += displacement
      let coldEnergy = coldCalculator.energy
      XCTAssertEqual(warmEnergy, coldEnergy, accuracy: 0.01)
    }
  }
  
  func testDifferentCalculator() throws {
    xTB_Environment.verbosity = .muted
    let calculator1 = SeedTests.createCalculator()
    let calculator2 = SeedTests.createCalculator()
    
    _ = calculator1.energy
    let seed = try XCTUnwrap(calculator1.seed)
    
    var bodyExecuted = false
    XCTAssertThrowsError(
      try calculator2.withSeed(seed) {
        bodyExecuted = true
      }
    ) { error in
      XCTAssertEqual(error as? xTB_Seed.Error, .differentCalculator)
    }
    XCTAssertFalse(bodyExecuted)
  }
  
  func testUnusedSeed() throws {
    xTB_Environment.verbosity = .muted
    let calculator = SeedTests.createCalculator()
    
    _ = calculator.energy
    let seed = try XCTUnwrap(calculator.seed)
    
    // The energy is already cached, so no singlepoint is computed.
    XCTAssertThrowsError(
      try calculator.withSeed(seed) {
        _ = calculator.energy
      }
    ) { error in
      XCTAssertEqual(error as? xTB_Seed.Error, .unused)
    }
    
    // The seed was forgotten, and does not leak into the next singlepoint.
    calculator.molecule.positions[0].x += 0.001
    XCTAssertNoThrow(
      try calculator.withSeed(seed) {
        _ = calculator.energy
      })
  }
  
  // When the test was written:
  // - warm: [7, 3, 7] SCC iterations
  // - cold: [12, 13, 12] SCC iterations
  
  func testWarmStartIterations() throws {
    xTB_Environment.verbosity = .muted
    let warmCalculator = SeedTests.createCalculator()
    let coldCalculator = SeedTests.createCalculator()
    
    _ = warmCalculator.energy
    let seed = try XCTUnwrap(warmCalculator.seed)
    
    let originalPositions = warmCalculator.molecule.positions
    let displacements: [Float] = [-0.001, 0.000, 0.001]
    for displacement in displacements {
      var positions = originalPositions
      positions[0].x += displacement

      // Warm start from the seed.
      let warmIterations = try warmCalculator.withSeed(seed) {
        warmCalculator.molecule.positions = positions
        return try SeedTests.countIterations {
          _ = warmCalculator.energy
        }
      }

      // Same calculator, outside the scope. Assigning the positions again
      // invalidates the cached singlepoint.
      warmCalculator.molecule.positions = positions
      let coldIterations1 = try SeedTests.countIterations {
        _ = warmCalculator.energy
      }

      // A calculator that never had a seed.
      coldCalculator.molecule.positions = positions
      let coldIterations2 = try SeedTests.countIterations {
        _ = coldCalculator.energy
      }
      
      XCTAssertGreaterThanOrEqual(warmIterations, 1)
      XCTAssertLessThan(warmIterations, 20)
      
      XCTAssertGreaterThanOrEqual(coldIterations1, 5)
      XCTAssertLessThan(coldIterations1, 40)
      
      print("iterations: (\(warmIterations), \(coldIterations1))")
      
      XCTAssertLessThan(warmIterations, coldIterations1)
      XCTAssertLessThan(warmIterations, coldIterations2)
      XCTAssertEqual(coldIterations1, coldIterations2)
    }
  }

  func testWarmStartDoesNotChain() throws {
    xTB_Environment.verbosity = .muted
    let calculator = SeedTests.createCalculator()
    _ = calculator.energy
    let seed = try XCTUnwrap(calculator.seed)
    
    let originalPositions = calculator.molecule.positions
    let displacements: [Float] = [-0.001, 0.000, 0.001]
    func countIterations(displacement: Float) throws -> Int {
      calculator.molecule.positions = originalPositions
      calculator.molecule.positions[0].x += displacement
      return try SeedTests.countIterations {
        _ = calculator.energy
      }
    }
    
    // Several singlepoints inside one scope.
    let sharedScopeIterations = try calculator.withSeed(seed) {
      try displacements.map(countIterations(displacement:))
    }
    
    // Each singlepoint inside its own scope.
    let separateScopeIterations = try displacements.map { displacement in
      try calculator.withSeed(seed) {
        try countIterations(displacement: displacement)
      }
    }
    
    // Every singlepoint starts from an identical copy of the seed, so the
    // counts match exactly. If the singlepoints chained, the undisplaced
    // geometry would start from the +0.001 nm solution instead of the seed.
    XCTAssertEqual(sharedScopeIterations, separateScopeIterations)
    
    // If these cases have the same iteration count, the test cannot function.
    XCTAssertNotEqual(sharedScopeIterations[0], sharedScopeIterations[1])
  }
}

/// Redirects file descriptor 1 of the C runtime used by libgfortran, which is
/// where xTB's console output goes.
///
/// On Windows, the MSYS2 MINGW64 build of libgfortran writes through
/// msvcrt.dll. Swift uses ucrtbase.dll, which has a separate file descriptor
/// table. So the redirection must call into msvcrt.dll. A build from the
/// UCRT64 environment would share Swift's runtime instead.
struct StdoutRedirection {
  #if os(Windows)
  typealias OpenFunction = @convention(c) (
    UnsafePointer<CChar>?, Int32, Int32) -> Int32
  typealias DupFunction = @convention(c) (Int32) -> Int32
  typealias Dup2Function = @convention(c) (Int32, Int32) -> Int32
  typealias CloseFunction = @convention(c) (Int32) -> Int32

  static func loadSymbol<T>(name: String, type: T.Type = T.self) -> T {
    let module = LoadLibraryA("msvcrt.dll")
    guard let symbol = GetProcAddress(module, name) else {
      fatalError("Could not load \(name) from msvcrt.dll.")
    }
    return unsafeBitCast(symbol, to: type)
  }

  let close: CloseFunction = loadSymbol(name: "_close")
  let dup2: Dup2Function = loadSymbol(name: "_dup2")
  #endif

  let fileDescriptor: Int32
  let savedStdout: Int32

  init(path: String) {
    #if os(Windows)
    let open: OpenFunction = Self.loadSymbol(name: "_open")
    let dup: DupFunction = Self.loadSymbol(name: "_dup")
    let _O_WRONLY: Int32 = 0x0001
    let _O_TRUNC: Int32 = 0x0200
    fileDescriptor = open(path, _O_WRONLY | _O_TRUNC, 0)
    savedStdout = dup(1)
    #else
    fileDescriptor = open(path, O_WRONLY | O_TRUNC)
    savedStdout = dup(1)
    fflush(stdout)
    #endif
    guard fileDescriptor >= 0, savedStdout >= 0 else {
      fatalError("Could not redirect stdout.")
    }
    _ = dup2(fileDescriptor, 1)
  }

  func restore() {
    #if !os(Windows)
    fflush(stdout)
    #endif
    _ = dup2(savedStdout, 1)
    _ = close(savedStdout)
    _ = close(fileDescriptor)
  }
}
