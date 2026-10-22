import XCTest
import xTB

final class PrecisionTests: XCTestCase {
  static func createCalculator(
    hamiltonian: xTB_Hamiltonian = .tightBinding
  ) -> xTB_Calculator {
    let system: [SIMD4<Float>] = diamondSystem122
    var calculatorDesc = xTB_CalculatorDescriptor()
    calculatorDesc.atomicNumbers = system.map { UInt8($0.w) }
    calculatorDesc.positions = system.map {
      SIMD3<Double>(SIMD3($0.x, $0.y, $0.z))
    }
    calculatorDesc.hamiltonian = hamiltonian
    return xTB_Calculator(descriptor: calculatorDesc)
  }

  func testDefaultPrecision() throws {
    xTB_Environment.verbosity = .muted
    let calculator = PrecisionTests.createCalculator()
    XCTAssertEqual(calculator.precision, .fp64)

    guard xTB_Library.supportsFP32 else {
      // Libraries without the SCC precision API report nothing.
      XCTAssertNil(calculator.precisionReport)
      return
    }
    let report = try XCTUnwrap(calculator.precisionReport)
    XCTAssertEqual(report.precision, .fp64)
    XCTAssertEqual(report.fp32Iterations, 0)
    XCTAssertFalse(report.restartedInFP64)
    XCTAssertGreaterThan(report.iterations, 0)
  }

  func testFP32() throws {
    guard xTB_Library.supportsFP32 else {
      throw XCTSkip("The xTB library has no SCC precision API.")
    }
    xTB_Environment.verbosity = .muted
    let fp64Calculator = PrecisionTests.createCalculator()
    let fp32Calculator = PrecisionTests.createCalculator()
    fp32Calculator.precision = .fp32()

    // Within the production accuracy gate, 5e-7 Hartree (2.2e-3 zJ; energies
    // here are in zJ). The worst case measured over the test set is ~2e-7.
    XCTAssertEqual(fp32Calculator.energy, fp64Calculator.energy, accuracy: 2e-3)

    let report = try XCTUnwrap(fp32Calculator.precisionReport)
    XCTAssertEqual(report.precision, .fp32(
      switchThreshold: 1e-3, maximumIterations: 160, refined: false))
    XCTAssertGreaterThan(report.fp32Iterations, 0)
    XCTAssertFalse(report.restartedInFP64)
  }

  func testRefinedFP32() throws {
    guard xTB_Library.supportsFP32 else {
      throw XCTSkip("The xTB library has no SCC precision API.")
    }
    xTB_Environment.verbosity = .muted
    let fp64Calculator = PrecisionTests.createCalculator()
    let refinedCalculator = PrecisionTests.createCalculator()
    refinedCalculator.precision = .fp32(refined: true)

    // About 1e-12 Hartree from FP64; allow 1e-9 Hartree (4.4e-6 zJ).
    XCTAssertEqual(
      refinedCalculator.energy, fp64Calculator.energy, accuracy: 4.4e-6)
    let report = try XCTUnwrap(refinedCalculator.precisionReport)
    XCTAssertEqual(report.precision, .fp32(refined: true))
  }

  func testChangingPrecision() throws {
    guard xTB_Library.supportsFP32 else {
      throw XCTSkip("The xTB library has no SCC precision API.")
    }
    xTB_Environment.verbosity = .muted
    let calculator = PrecisionTests.createCalculator()
    let fp64Energy = calculator.energy
    XCTAssertNotNil(calculator.seed)

    // Changing the precision invalidates the singlepoint.
    calculator.precision = .fp32(switchThreshold: 1e-2)
    XCTAssertNil(calculator.seed)
    _ = calculator.energy
    XCTAssertEqual(
      calculator.precisionReport?.precision, .fp32(switchThreshold: 1e-2))

    // Back to FP64: the same result as before. Bit for bit with one OpenMP
    // thread; with several (test.sh), xTB's threaded reductions differ from
    // run to run at about 1e-15 relative.
    calculator.precision = .fp64
    XCTAssertEqual(calculator.energy, fp64Energy, accuracy: 1e-6)
    XCTAssertEqual(calculator.precisionReport?.precision, .fp64)
  }

  func testForceField() throws {
    xTB_Environment.verbosity = .muted
    let calculator = PrecisionTests.createCalculator(hamiltonian: .forceField)
    XCTAssertEqual(calculator.precision, .fp64)
    XCTAssertNil(calculator.precisionReport)
  }
}
