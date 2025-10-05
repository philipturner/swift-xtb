import XCTest
import xTB

final class GFNFF_Tests: XCTestCase {
  func testDiamondSystem122() throws {
    // Set the environment verbosity.
    xTB_Environment.verbosity = .muted

    // Select the system.
    let system: [SIMD4<Float>] = diamondSystem122

    // Create the calculator.
    var calculatorDesc = xTB_CalculatorDescriptor()
    calculatorDesc.atomicNumbers = system.map { UInt8($0.w) }
    calculatorDesc.positions = system.map {
      SIMD3($0.x, $0.y, $0.z)
    }
    calculatorDesc.hamiltonian = .forceField
    let calculator = xTB_Calculator(descriptor: calculatorDesc)
    XCTAssertEqual(calculator.molecule.atomicNumbers.count, 82)

    let energy = calculator.energy
    let formattedEnergy = String(format: "%.3f", energy)
    print("energy:", formattedEnergy, "zJ")
  }
}
