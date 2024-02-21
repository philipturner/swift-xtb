import XCTest
import Mechanosynthesis
import Numerics

final class AnsatzTests: XCTestCase {
  static func checkFragments(
    _ waveFunction: WaveFunction,
    _ expectedCount: Int
  ) {
    print("fragment count:", waveFunction.cellValues.count)
    XCTAssertGreaterThanOrEqual(waveFunction.cellValues.count, expectedCount)
    XCTAssertLessThan(waveFunction.cellValues.count, expectedCount * 2)
  }
  
  static func queryRadius(
    waveFunction: WaveFunction,
    nucleusPosition: SIMD3<Float> = .zero
  ) -> Float {
    var sum: Double = .zero
    let octree = waveFunction.octree
    for nodeID in octree.linkedList.indices {
      if octree.linkedList[nodeID].childCount > 0 {
        continue
      }
      
      let metadata = octree.metadata[nodeID]
      var x = SIMD8<Float>(0, 1, 0, 1, 0, 1, 0, 1) * 0.5 - 0.25
      var y = SIMD8<Float>(0, 0, 1, 1, 0, 0, 1, 1) * 0.5 - 0.25
      var z = SIMD8<Float>(0, 0, 0, 0, 1, 1, 1, 1) * 0.5 - 0.25
      x = x * metadata.w + metadata.x
      y = y * metadata.w + metadata.y
      z = z * metadata.w + metadata.z
      x -= nucleusPosition.x
      y -= nucleusPosition.y
      z -= nucleusPosition.z
      
      let Ψ = waveFunction.cellValues[nodeID]
      let r = (x * x + y * y + z * z).squareRoot()
      let ΨrΨ = Ψ * r * Ψ
      
      let d3r = metadata.w * metadata.w * metadata.w
      sum += Double(ΨrΨ.sum() / 8 * d3r)
    }
    
    let output = Float(sum) * 2 / 3
    print("expectation radius:", output)
    return output
  }
  
  func testHydrogen() throws {
    print()
    print("testHydrogen")
    
    var descriptor = AnsatzDescriptor()
    descriptor.atomicNumbers = [1]
    descriptor.fragmentCount = 1000
    descriptor.positions = [.zero]
    descriptor.sizeExponent = 4
    
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
    Self.checkFragments(hydrogenUp.spinUpWaveFunctions[0], 1000)
    XCTAssertEqual(1, Self.queryRadius(
      waveFunction: hydrogenUp.spinUpWaveFunctions[0]), accuracy: 0.02)
    
    // Test a neutral hydrogen atom (negative spin).
    descriptor.netCharges = [0]
    descriptor.netSpinPolarizations = [-1]
    let hydrogenDown = Ansatz(descriptor: descriptor)
    XCTAssertEqual(hydrogenDown.spinDownWaveFunctions.count, 1)
    XCTAssertEqual(hydrogenDown.spinNeutralWaveFunctions.count, 0)
    XCTAssertEqual(hydrogenDown.spinUpWaveFunctions.count, 0)
    Self.checkFragments(hydrogenDown.spinDownWaveFunctions[0], 1000)
    XCTAssertEqual(1, Self.queryRadius(
      waveFunction: hydrogenDown.spinDownWaveFunctions[0]), accuracy: 0.02)
    
    // Test a hydride anion (singlet state).
    descriptor.netCharges = [-1]
    descriptor.netSpinPolarizations = [0]
    let hydrideSinglet = Ansatz(descriptor: descriptor)
    XCTAssertEqual(hydrideSinglet.spinDownWaveFunctions.count, 0)
    XCTAssertEqual(hydrideSinglet.spinNeutralWaveFunctions.count, 1)
    XCTAssertEqual(hydrideSinglet.spinUpWaveFunctions.count, 0)
    Self.checkFragments(hydrideSinglet.spinNeutralWaveFunctions[0], 1000)
    XCTAssertEqual(1, Self.queryRadius(
      waveFunction: hydrideSinglet.spinNeutralWaveFunctions[0]), accuracy: 0.02)
    
    // Test a hydride anion (triplet state).
    descriptor.atomicNumbers = [1]
    descriptor.netCharges = [-1]
    descriptor.netSpinPolarizations = [2]
    let hydrideTriplet = Ansatz(descriptor: descriptor)
    XCTAssertEqual(hydrideTriplet.spinDownWaveFunctions.count, 0)
    XCTAssertEqual(hydrideTriplet.spinNeutralWaveFunctions.count, 0)
    XCTAssertEqual(hydrideTriplet.spinUpWaveFunctions.count, 2)
    Self.checkFragments(hydrideTriplet.spinUpWaveFunctions[0], 1000)
    Self.checkFragments(hydrideTriplet.spinUpWaveFunctions[1], 1000)
    XCTAssertEqual(1, Self.queryRadius(
      waveFunction: hydrideTriplet.spinUpWaveFunctions[0]), accuracy: 0.02)
    XCTAssertEqual(2, Self.queryRadius(
      waveFunction: hydrideTriplet.spinUpWaveFunctions[1]), accuracy: 0.08)
  }
  
  func testLithium() throws {
    print()
    print("testLithium")
    
    var descriptor = AnsatzDescriptor()
    descriptor.atomicNumbers = [3]
    descriptor.fragmentCount = 1000
    descriptor.positions = [.zero]
    descriptor.sizeExponent = 4
    
    // Test an ion.
    descriptor.netCharges = [+1]
    descriptor.netSpinPolarizations = [0]
    let lithiumIon = Ansatz(descriptor: descriptor)
    XCTAssertEqual(lithiumIon.spinDownWaveFunctions.count, 0)
    XCTAssertEqual(lithiumIon.spinNeutralWaveFunctions.count, 1)
    XCTAssertEqual(lithiumIon.spinUpWaveFunctions.count, 0)
    Self.checkFragments(lithiumIon.spinNeutralWaveFunctions[0], 1000)
    XCTAssertEqual(0.414, Self.queryRadius(
      waveFunction: lithiumIon.spinNeutralWaveFunctions[0]), accuracy: 0.02)
    
    // Test a neutral atom.
    descriptor.netCharges = [0]
    descriptor.netSpinPolarizations = [-1]
    let lithiumNeutral = Ansatz(descriptor: descriptor)
    XCTAssertEqual(lithiumNeutral.spinDownWaveFunctions.count, 1)
    XCTAssertEqual(lithiumNeutral.spinNeutralWaveFunctions.count, 1)
    XCTAssertEqual(lithiumNeutral.spinUpWaveFunctions.count, 0)
    Self.checkFragments(lithiumNeutral.spinDownWaveFunctions[0], 1000)
    Self.checkFragments(lithiumNeutral.spinNeutralWaveFunctions[0], 1000)
    XCTAssertEqual(0.414, Self.queryRadius(
      waveFunction: lithiumNeutral.spinNeutralWaveFunctions[0]), accuracy: 0.02)
    XCTAssertEqual(2, Self.queryRadius(
      waveFunction: lithiumNeutral.spinDownWaveFunctions[0]), accuracy: 0.08)
    
    // Test a spin-3/2 atom.
    descriptor.netCharges = [0]
    descriptor.netSpinPolarizations = [3]
    let lithiumPolarized = Ansatz(descriptor: descriptor)
    XCTAssertEqual(lithiumPolarized.spinDownWaveFunctions.count, 0)
    XCTAssertEqual(lithiumPolarized.spinNeutralWaveFunctions.count, 0)
    XCTAssertEqual(lithiumPolarized.spinUpWaveFunctions.count, 3)
    Self.checkFragments(lithiumPolarized.spinUpWaveFunctions[0], 1000)
    Self.checkFragments(lithiumPolarized.spinUpWaveFunctions[1], 1000)
    Self.checkFragments(lithiumPolarized.spinUpWaveFunctions[2], 1000)
    XCTAssertEqual(0.333, Self.queryRadius(
      waveFunction: lithiumPolarized.spinUpWaveFunctions[0]), accuracy: 0.02)
    XCTAssertEqual(1.414, Self.queryRadius(
      waveFunction: lithiumPolarized.spinUpWaveFunctions[1]), accuracy: 0.03)
    XCTAssertEqual(1.414 * 5 / 6, Self.queryRadius(
      waveFunction: lithiumPolarized.spinUpWaveFunctions[2]), accuracy: 0.03)
  }
  
  // Test the nitrogen dimer.
  func testNitrogen() throws {
    print()
    print("testNitrogen")
    
    var descriptor = AnsatzDescriptor()
    descriptor.atomicNumbers = [7, 7]
    descriptor.fragmentCount = 1000
    descriptor.positions = [SIMD3(-1, 0, 0), SIMD3(1, 0, 0)]
    descriptor.sizeExponent = 4
    
    // Test N+ and N-.
    descriptor.netCharges = [+1, -1]
    descriptor.netSpinPolarizations = [0, 0]
    let nitrogenCharged = Ansatz(descriptor: descriptor)
    XCTAssertEqual(nitrogenCharged.spinDownWaveFunctions.count, 0)
    XCTAssertEqual(nitrogenCharged.spinNeutralWaveFunctions.count, 7)
    XCTAssertEqual(nitrogenCharged.spinUpWaveFunctions.count, 0)
    for i in 0..<7 {
      Self.checkFragments(nitrogenCharged.spinNeutralWaveFunctions[i], 1000)
    }
    XCTAssertEqual(1 / 6.414, Self.queryRadius(
      waveFunction: nitrogenCharged.spinNeutralWaveFunctions[0],
      nucleusPosition: SIMD3(-1, 0, 0)), accuracy: 5e-3)
    XCTAssertEqual(0.667, Self.queryRadius(
      waveFunction: nitrogenCharged.spinNeutralWaveFunctions[1],
      nucleusPosition: SIMD3(-1, 0, 0)), accuracy: 0.08)
    XCTAssertEqual(0.667 * 5 / 6, Self.queryRadius(
      waveFunction: nitrogenCharged.spinNeutralWaveFunctions[2],
      nucleusPosition: SIMD3(-1, 0, 0)), accuracy: 0.08)
    XCTAssertEqual(1 / 6.414, Self.queryRadius(
      waveFunction: nitrogenCharged.spinNeutralWaveFunctions[3],
      nucleusPosition: SIMD3(1, 0, 0)), accuracy: 5e-3)
    XCTAssertEqual(2 / 1.449, Self.queryRadius(
      waveFunction: nitrogenCharged.spinNeutralWaveFunctions[4],
      nucleusPosition: SIMD3(1, 0, 0)), accuracy: 0.04)
    XCTAssertEqual(2 / 1.449 * 5 / 6, Self.queryRadius(
      waveFunction: nitrogenCharged.spinNeutralWaveFunctions[5],
      nucleusPosition: SIMD3(1, 0, 0)), accuracy: 0.04)
    XCTAssertEqual(2 / 1.449 * 5 / 6, Self.queryRadius(
      waveFunction: nitrogenCharged.spinNeutralWaveFunctions[6],
      nucleusPosition: SIMD3(1, 0, 0)), accuracy: 0.04)
    
    // Test N(3+) anions.
    // - Literature value for ion atomic radius is 140 pm, or 2.65 Bohr.
    // - The ansatz produces 2.00 Bohr - not bad.
    descriptor.netCharges = [-3, -3]
    descriptor.netSpinPolarizations = [0, 0]
    let nitride = Ansatz(descriptor: descriptor)
    XCTAssertEqual(nitride.spinDownWaveFunctions.count, 0)
    XCTAssertEqual(nitride.spinNeutralWaveFunctions.count, 10)
    XCTAssertEqual(nitride.spinUpWaveFunctions.count, 0)
    for i in 0..<10 {
      Self.checkFragments(nitride.spinNeutralWaveFunctions[i], 1000)
    }
    XCTAssertEqual(1 / 6.414, Self.queryRadius(
      waveFunction: nitride.spinNeutralWaveFunctions[0],
      nucleusPosition: SIMD3(-1, 0, 0)), accuracy: 5e-3)
    XCTAssertEqual(2, Self.queryRadius(
      waveFunction: nitride.spinNeutralWaveFunctions[1],
      nucleusPosition: SIMD3(-1, 0, 0)), accuracy: 0.06)
    XCTAssertEqual(2 * 5 / 6, Self.queryRadius(
      waveFunction: nitride.spinNeutralWaveFunctions[2],
      nucleusPosition: SIMD3(-1, 0, 0)), accuracy: 0.06)
    XCTAssertEqual(2 * 5 / 6, Self.queryRadius(
      waveFunction: nitride.spinNeutralWaveFunctions[3],
      nucleusPosition: SIMD3(-1, 0, 0)), accuracy: 0.06)
    XCTAssertEqual(2 * 5 / 6, Self.queryRadius(
      waveFunction: nitride.spinNeutralWaveFunctions[4],
      nucleusPosition: SIMD3(-1, 0, 0)), accuracy: 0.06)
    XCTAssertEqual(1 / 6.414, Self.queryRadius(
      waveFunction: nitride.spinNeutralWaveFunctions[5],
      nucleusPosition: SIMD3(1, 0, 0)), accuracy: 5e-3)
    XCTAssertEqual(2, Self.queryRadius(
      waveFunction: nitride.spinNeutralWaveFunctions[6],
      nucleusPosition: SIMD3(1, 0, 0)), accuracy: 0.06)
    XCTAssertEqual(2 * 5 / 6, Self.queryRadius(
      waveFunction: nitride.spinNeutralWaveFunctions[7],
      nucleusPosition: SIMD3(1, 0, 0)), accuracy: 0.06)
    XCTAssertEqual(2 * 5 / 6, Self.queryRadius(
      waveFunction: nitride.spinNeutralWaveFunctions[8],
      nucleusPosition: SIMD3(1, 0, 0)), accuracy: 0.06)
    XCTAssertEqual(2 * 5 / 6, Self.queryRadius(
      waveFunction: nitride.spinNeutralWaveFunctions[9],
      nucleusPosition: SIMD3(1, 0, 0)), accuracy: 0.06)
    
    // Test the triplet state.
    descriptor.netCharges = [0, 0]
    descriptor.netSpinPolarizations = [1, 1]
    let nitrogenTriplet = Ansatz(descriptor: descriptor)
    XCTAssertEqual(nitrogenTriplet.spinDownWaveFunctions.count, 0)
    XCTAssertEqual(nitrogenTriplet.spinNeutralWaveFunctions.count, 6)
    XCTAssertEqual(nitrogenTriplet.spinUpWaveFunctions.count, 2)
    for i in 0..<6 {
      Self.checkFragments(nitrogenTriplet.spinNeutralWaveFunctions[i], 1000)
    }
    Self.checkFragments(nitrogenTriplet.spinUpWaveFunctions[0], 1000)
    Self.checkFragments(nitrogenTriplet.spinUpWaveFunctions[1], 1000)
    
    // Test an ROHF simulation that forms 3 covalent bonds.
    descriptor.netCharges = [0, 0]
    descriptor.netSpinPolarizations = [3, -3]
    let nitrogenROHF = Ansatz(descriptor: descriptor)
    XCTAssertEqual(nitrogenROHF.spinDownWaveFunctions.count, 3)
    XCTAssertEqual(nitrogenROHF.spinNeutralWaveFunctions.count, 4)
    XCTAssertEqual(nitrogenROHF.spinUpWaveFunctions.count, 3)
    for i in 0..<4 {
      Self.checkFragments(nitrogenROHF.spinNeutralWaveFunctions[i], 1000)
    }
    for i in 0..<3 {
      Self.checkFragments(nitrogenROHF.spinDownWaveFunctions[i], 1000)
      Self.checkFragments(nitrogenROHF.spinUpWaveFunctions[i], 1000)
    }
  }
  
  // Test chromium.
  // - Test the limitation of GFN-xTB, which doesn't support the true spin
  //   polarization of chromium d-orbitals. Instead, it does 3d^4, 4s^2.
  func testChromium() throws {
    print()
    print("testChromium")
    
    var descriptor = AnsatzDescriptor()
    descriptor.atomicNumbers = [24]
    descriptor.fragmentCount = 1000
    descriptor.positions = [.zero]
    descriptor.sizeExponent = 4
    
    // Test the 3+ oxidation state, with minimal spin polarization.
    descriptor.netCharges = [+3]
    descriptor.netSpinPolarizations = [1]
    let chromiumIon = Ansatz(descriptor: descriptor)
    XCTAssertEqual(chromiumIon.spinDownWaveFunctions.count, 0)
    XCTAssertEqual(chromiumIon.spinNeutralWaveFunctions.count, 10)
    XCTAssertEqual(chromiumIon.spinUpWaveFunctions.count, 1)
    for i in 0..<10 {
      Self.checkFragments(chromiumIon.spinNeutralWaveFunctions[i], 1000)
    }
    Self.checkFragments(chromiumIon.spinUpWaveFunctions[0], 1000)
    XCTAssertEqual(1 / 23.414, Self.queryRadius(
      waveFunction: chromiumIon.spinNeutralWaveFunctions[0]), accuracy: 2e-3)
    XCTAssertEqual(2 / 16.828, Self.queryRadius(
      waveFunction: chromiumIon.spinNeutralWaveFunctions[1]), accuracy: 4e-3)
    XCTAssertEqual(2 / 16.828 * 5 / 6, Self.queryRadius(
      waveFunction: chromiumIon.spinNeutralWaveFunctions[2]), accuracy: 4e-3)
    XCTAssertEqual(2 / 16.828 * 5 / 6, Self.queryRadius(
      waveFunction: chromiumIon.spinNeutralWaveFunctions[3]), accuracy: 4e-3)
    XCTAssertEqual(2 / 16.828 * 5 / 6, Self.queryRadius(
      waveFunction: chromiumIon.spinNeutralWaveFunctions[4]), accuracy: 4e-3)
    XCTAssertEqual(0.375, Self.queryRadius(
      waveFunction: chromiumIon.spinNeutralWaveFunctions[5]), accuracy: 0.01)
    XCTAssertEqual(0.375 * 12 / 13, Self.queryRadius(
      waveFunction: chromiumIon.spinNeutralWaveFunctions[6]), accuracy: 0.01)
    XCTAssertEqual(0.375 * 12 / 13, Self.queryRadius(
      waveFunction: chromiumIon.spinNeutralWaveFunctions[7]), accuracy: 0.01)
    XCTAssertEqual(0.375 * 12 / 13, Self.queryRadius(
      waveFunction: chromiumIon.spinNeutralWaveFunctions[8]), accuracy: 0.01)
    XCTAssertEqual(0.2906231, Self.queryRadius(
      waveFunction: chromiumIon.spinUpWaveFunctions[0]), accuracy: 0.01)
    XCTAssertEqual(4 / 4.414, Self.queryRadius(
      waveFunction: chromiumIon.spinNeutralWaveFunctions[9]), accuracy: 0.02)
    
    // Test an Aufbau polarized atom.
    descriptor.netCharges = [0]
    descriptor.netSpinPolarizations = [4]
    let chromiumAufbau = Ansatz(descriptor: descriptor)
    XCTAssertEqual(chromiumAufbau.spinDownWaveFunctions.count, 0)
    XCTAssertEqual(chromiumAufbau.spinNeutralWaveFunctions.count, 10)
    XCTAssertEqual(chromiumAufbau.spinUpWaveFunctions.count, 4)
    for i in 0..<10 {
      Self.checkFragments(chromiumAufbau.spinNeutralWaveFunctions[i], 1000)
    }
    for i in 0..<4 {
      Self.checkFragments(chromiumAufbau.spinUpWaveFunctions[i], 1000)
    }
    
    // Test an actual polarized atom.
    descriptor.netCharges = [0]
    descriptor.netSpinPolarizations = [6]
    let chromiumActual = Ansatz(descriptor: descriptor)
    XCTAssertEqual(chromiumActual.spinDownWaveFunctions.count, 0)
    XCTAssertEqual(chromiumActual.spinNeutralWaveFunctions.count, 9)
    XCTAssertEqual(chromiumActual.spinUpWaveFunctions.count, 6)
    for i in 0..<9 {
      Self.checkFragments(chromiumActual.spinNeutralWaveFunctions[i], 1000)
    }
    for i in 0..<6 {
      Self.checkFragments(chromiumActual.spinUpWaveFunctions[i], 1000)
    }
  }
  
  // Test group (IV) atoms.
  // - Compare to graphs and spreadsheet data from 'AnsatzExperiment'.
  
  // Then, optimize the code for generating the initial guess. Avoid use of
  // multiple cores and AMX for the foreseeable future, to get the
  // highest-quality data about compute cost.
}
