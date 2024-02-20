import XCTest
import Mechanosynthesis
import Numerics

final class AnsatzTests: XCTestCase {
  func testHydrogen() throws {
    var descriptor = AnsatzDescriptor()
    descriptor.atomicNumbers = [1]
    descriptor.fragmentCount = 1000
    descriptor.positions = [.zero]
    descriptor.sizeExponent = 4
    
    func checkFragmentCount(_ waveFunction: WaveFunction) {
      XCTAssertGreaterThanOrEqual(waveFunction.cellValues.count, 1000)
      XCTAssertLessThan(waveFunction.cellValues.count, 2000)
    }
    
    // Test a proton.
    descriptor.netCharges = [+1]
    descriptor.netSpinPolarizations = [0]
    let proton = Ansatz(descriptor: descriptor)
    XCTAssertEqual(proton.spinDownWaveFunctions.count, 0)
    XCTAssertEqual(proton.spinNeutralWaveFunctions.count, 0)
    XCTAssertEqual(proton.spinUpWaveFunctions.count, 0)
    
    // Test a neutral hydrogen atom (positive spin).
    descriptor.netCharges = [0]
    descriptor.netSpinPolarizations = [1]
    let hydrogenUp = Ansatz(descriptor: descriptor)
    XCTAssertEqual(hydrogenUp.spinDownWaveFunctions.count, 0)
    XCTAssertEqual(hydrogenUp.spinNeutralWaveFunctions.count, 0)
    XCTAssertEqual(hydrogenUp.spinUpWaveFunctions.count, 1)
    checkFragmentCount(hydrogenUp.spinUpWaveFunctions[0])
    
    // Test a neutral hydrogen atom (negative spin).
    descriptor.netCharges = [0]
    descriptor.netSpinPolarizations = [-1]
    let hydrogenDown = Ansatz(descriptor: descriptor)
    XCTAssertEqual(hydrogenDown.spinDownWaveFunctions.count, 1)
    XCTAssertEqual(hydrogenDown.spinNeutralWaveFunctions.count, 0)
    XCTAssertEqual(hydrogenDown.spinUpWaveFunctions.count, 0)
    checkFragmentCount(hydrogenDown.spinDownWaveFunctions[0])
    
    // Test a hydride anion.
    descriptor.netCharges = [-1]
    descriptor.netSpinPolarizations = [0]
    let hydride = Ansatz(descriptor: descriptor)
    XCTAssertEqual(hydride.spinDownWaveFunctions.count, 0)
    XCTAssertEqual(hydride.spinNeutralWaveFunctions.count, 1)
    XCTAssertEqual(hydride.spinUpWaveFunctions.count, 0)
    checkFragmentCount(hydride.spinNeutralWaveFunctions[0])
    
    // TODO: Test expectation values now.
  }
  
  func testHelium() throws {
    
  }
}
