import Dispatch
import Foundation
#if os(Windows)
import WinSDK
#endif
import XCTest
@testable import xTB

// `xTB_Environment.threadCount`. The OpenMP runtime reads OMP_NUM_THREADS when
// it loads; the count is per thread afterward, and swift-xtb sets it on the
// thread that starts each calculation. The team size of a parallel region is
// measured directly: an empty region started through libgomp's GOMP_parallel
// reports omp_get_num_threads. A thread that never ran a calculation shows the
// runtime's own setting, whatever other tests set.
final class OpenMPTests: XCTestCase {
  // MARK: - Utilities

  /// A symbol of the loaded OpenMP runtime (libgomp), or nil.
  static func openMPSymbol(_ name: String) -> UnsafeMutableRawPointer? {
    #if os(Windows)
    guard let module = GetModuleHandleA("libgomp-1.dll"),
          let symbol = GetProcAddress(module, name) else {
      return nil
    }
    return unsafeBitCast(symbol, to: UnsafeMutableRawPointer?.self)
    #else
    return dlsym(UnsafeMutableRawPointer(bitPattern: -2), name) // RTLD_DEFAULT
    #endif
  }

  nonisolated(unsafe) static var lastTeamSize: Int32 = -1

  /// The number of threads in a parallel region started on the calling
  /// thread (libgomp only).
  static func teamSize() throws -> Int {
    typealias Region = @convention(c) (UnsafeMutableRawPointer?) -> Void
    typealias Parallel = @convention(c) (
      Region, UnsafeMutableRawPointer?, UInt32, UInt32) -> Void
    let parallel = try XCTUnwrap(
      openMPSymbol("GOMP_parallel"), "GOMP_parallel not found (not libgomp?)")
    lastTeamSize = -1
    unsafeBitCast(parallel, to: Parallel.self)({ _ in
      typealias Count = @convention(c) () -> Int32
      let count = unsafeBitCast(
        OpenMPTests.openMPSymbol("omp_get_num_threads")!, to: Count.self)
      OpenMPTests.lastTeamSize = count()
    }, nil, 0, 0)
    return Int(lastTeamSize)
  }

  /// The outcome of a closure run on another thread, handed back to the
  /// waiting one (the closures that run there must be @Sendable, and so can't
  /// write a captured variable).
  final class Outcome<T>: @unchecked Sendable {
    private let lock = NSLock()
    private var result: Result<T, Error>?
    
    func set(_ body: @Sendable () throws -> T) {
      let result = Result { try body() }
      lock.withLock { self.result = result }
    }
    
    func get() throws -> T {
      try lock.withLock { result! }.get()
    }
  }
  
  /// Runs `body` on a new thread and waits for it.
  static func onNewThread<T>(
    _ body: @escaping @Sendable () throws -> T
  ) throws -> T {
    let outcome = Outcome<T>()
    let done = DispatchSemaphore(value: 0)
    let thread = Thread {
      outcome.set(body)
      done.signal()
    }
    thread.start()
    done.wait()
    return try outcome.get()
  }
  
  /// Runs `body` on a Dispatch worker and waits for it.
  static func onDispatchWorker<T>(
    _ body: @escaping @Sendable () throws -> T
  ) throws -> T {
    let outcome = Outcome<T>()
    let done = DispatchSemaphore(value: 0)
    DispatchQueue.global().async {
      outcome.set(body)
      done.signal()
    }
    done.wait()
    return try outcome.get()
  }
  
  static func createCalculator() -> xTB_Calculator {
    let system: [SIMD4<Float>] = diamondSystem122
    var calculatorDesc = xTB_CalculatorDescriptor()
    calculatorDesc.atomicNumbers = system.map { UInt8($0.w) }
    calculatorDesc.positions = system.map {
      SIMD3<Double>(SIMD3($0.x, $0.y, $0.z))
    }
    calculatorDesc.hamiltonian = .tightBinding
    return xTB_Calculator(descriptor: calculatorDesc)
  }

  /// The runtime's own setting: OMP_NUM_THREADS at launch, else every core.
  static var launchThreadCount: Int {
    if let value = ProcessInfo.processInfo.environment["OMP_NUM_THREADS"],
       let count = Int(value.split(separator: ",").first ?? "") {
      return count
    }
    return xTB_Environment.maximumThreadCount
  }

  // MARK: - Tests

  /// Without the API, the environment at launch decides.
  func testDefault() throws {
    xTB_Environment.verbosity = .muted
    let team = try OpenMPTests.onNewThread { try OpenMPTests.teamSize() }
    XCTAssertEqual(team, OpenMPTests.launchThreadCount)

    // The nominal core count is what the runtime uses without the variable.
    typealias Count = @convention(c) () -> Int32
    let processors = try XCTUnwrap(OpenMPTests.openMPSymbol("omp_get_num_procs"))
    XCTAssertEqual(
      xTB_Environment.maximumThreadCount,
      Int(unsafeBitCast(processors, to: Count.self)()))
  }

  /// With the API, every calculation runs with the count, on any thread.
  func testThreadCount() throws {
    xTB_Environment.verbosity = .muted
    defer { xTB_Environment.threadCount = nil }
    let initial = OpenMPTests.launchThreadCount
    let (a, b) = initial == 2 ? (3, 4) : initial == 3 ? (2, 4) : (2, 3)

    // Setting the count doesn't touch any thread yet.
    xTB_Environment.threadCount = a
    XCTAssertEqual(xTB_Environment.threadCount, a)
    XCTAssertEqual(
      try OpenMPTests.onNewThread { try OpenMPTests.teamSize() }, initial)

    // A calculation sets it on the thread that starts it: this one, a new
    // thread, a Dispatch worker.
    let energy = OpenMPTests.createCalculator().energy
    XCTAssert(energy.isFinite && energy < 0)
    XCTAssertEqual(try OpenMPTests.teamSize(), a)
    for runner in [OpenMPTests.onNewThread, OpenMPTests.onDispatchWorker]
      as [(@escaping @Sendable () throws -> (Double, Int)) throws -> (Double, Int)] {
      let (energy, team) = try runner {
        let energy = OpenMPTests.createCalculator().energy
        return (energy, try OpenMPTests.teamSize())
      }
      XCTAssert(energy.isFinite && energy < 0)
      XCTAssertEqual(team, a)
    }

    // A new count applies from the next calculation on.
    xTB_Environment.threadCount = b
    XCTAssertEqual(try OpenMPTests.teamSize(), a)
    _ = OpenMPTests.createCalculator().energy
    XCTAssertEqual(try OpenMPTests.teamSize(), b)
    let team = try OpenMPTests.onDispatchWorker {
      _ = OpenMPTests.createCalculator().energy
      return try OpenMPTests.teamSize()
    }
    XCTAssertEqual(team, b)

    // nil stops the overrides: a thread that never ran a calculation keeps
    // the runtime's setting, one that did keeps its last count.
    xTB_Environment.threadCount = nil
    XCTAssertNil(xTB_Environment.threadCount)
    _ = OpenMPTests.createCalculator().energy
    XCTAssertEqual(try OpenMPTests.teamSize(), b)
    let fresh = try OpenMPTests.onNewThread {
      _ = OpenMPTests.createCalculator().energy
      return try OpenMPTests.teamSize()
    }
    XCTAssertEqual(fresh, initial)
  }

  /// The same singlepoint with 1 and 4 threads: OpenMP only reorders sums.
  func testEnergyIndependentOfThreads() throws {
    xTB_Environment.verbosity = .muted
    defer { xTB_Environment.threadCount = nil }

    // The molecule refers to its calculator without owning it: keep both.
    xTB_Environment.threadCount = 1
    let calculator1 = OpenMPTests.createCalculator()
    let energy1 = calculator1.energy
    let forces1 = calculator1.molecule.forces
    XCTAssertEqual(try OpenMPTests.teamSize(), 1)
    xTB_Environment.threadCount = 4
    let calculator4 = OpenMPTests.createCalculator()
    let energy4 = calculator4.energy
    let forces4 = calculator4.molecule.forces
    XCTAssertEqual(try OpenMPTests.teamSize(), 4)

    XCTAssertEqual(energy1, energy4, accuracy: 1e-3) // zJ (2e-7 Eh)
    XCTAssertEqual(forces1.count, forces4.count)
    for (f1, f4) in zip(forces1, forces4) {
      let difference = f1 - f4
      let maximum = max(abs(difference.x), abs(difference.y), abs(difference.z))
      XCTAssertLessThan(maximum, 1) // pN
    }
  }
  
  /// Each step between `threadCount` and the runtime, on its own (the Windows
  /// failure documented in xTB_Environment+OpenMP.swift showed none of them
  /// separately): a direct omp_set_num_threads, the library's view of the
  /// runtime, applyThreadCount alone, then creating a calculator and running
  /// a singlepoint.
  func testRuntimeLookup() throws {
    xTB_Environment.verbosity = .muted
    defer { xTB_Environment.threadCount = nil }
    _ = OpenMPTests.createCalculator().energy
    let initial = OpenMPTests.launchThreadCount
    let (a, b) = initial == 2 ? (3, 4) : initial == 3 ? (2, 4) : (2, 3)
    
    // The runtime honors omp_set_num_threads called directly.
    typealias Set = @convention(c) (Int32) -> Void
    let setter = try XCTUnwrap(OpenMPTests.openMPSymbol("omp_set_num_threads"))
    let setterAddress = Int(bitPattern: setter)
    let direct = try OpenMPTests.onNewThread {
      let setter = UnsafeMutableRawPointer(bitPattern: setterAddress)!
      unsafeBitCast(setter, to: Set.self)(Int32(a))
      return try OpenMPTests.teamSize()
    }
    XCTAssertEqual(direct, a, "The runtime ignored omp_set_num_threads.")
    
    // The library finds the same runtime.
    XCTAssertTrue(xTB_Library.isLoaded)
    XCTAssertEqual(
      xTB_Library.loadOpenMPSymbol("omp_set_num_threads"), setter,
      "The library found a different omp_set_num_threads.")
    
    // applyThreadCount alone, on a fresh thread.
    xTB_Environment.threadCount = b
    typealias Count = @convention(c) () -> Int32
    let maxAddress = Int(bitPattern: try XCTUnwrap(
      OpenMPTests.openMPSymbol("omp_get_max_threads")))
    let (seen, maximum, applied) = try OpenMPTests.onNewThread {
      let seen = xTB_Environment.threadCount
      xTB_Environment.applyThreadCount()
      let maxThreads = UnsafeMutableRawPointer(bitPattern: maxAddress)!
      let maximum = unsafeBitCast(maxThreads, to: Count.self)()
      return (seen, Int(maximum), try OpenMPTests.teamSize())
    }
    XCTAssertEqual(seen, b, "threadCount reads differently on another thread.")
    XCTAssertEqual(maximum, b, "applyThreadCount didn't set the count.")
    XCTAssertEqual(applied, b, "applyThreadCount didn't set the count.")
    
    // Creating a calculator (parameter load), then a singlepoint.
    let (afterLoad, afterSinglepoint) = try OpenMPTests.onNewThread {
      let calculator = OpenMPTests.createCalculator()
      let afterLoad = try OpenMPTests.teamSize()
      _ = calculator.energy
      return (afterLoad, try OpenMPTests.teamSize())
    }
    XCTAssertEqual(afterLoad, b, "Creating a calculator didn't set the count.")
    XCTAssertEqual(afterSinglepoint, b, "A singlepoint didn't keep the count.")
  }
}
