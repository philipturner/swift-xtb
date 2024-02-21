import XCTest
import Mechanosynthesis
import Numerics

final class AnsatzTests: XCTestCase {
  static func checkFragments(
    _ waveFunction: WaveFunction,
    _ expectedCount: Int,
    _ upperExpectedCount: Int? = nil
  ) {
    XCTAssertGreaterThanOrEqual(waveFunction.cellValues.count, expectedCount)
    XCTAssertLessThan(
      waveFunction.cellValues.count, upperExpectedCount ?? expectedCount * 2)
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
    
    return Float(sum) * 2 / 3
  }
  
  func testHydrogen() throws {
    var descriptor = AnsatzDescriptor()
    descriptor.atomicNumbers = [1]
    descriptor.fragmentCount = 1000
    descriptor.positions = [.zero]
    descriptor.sizeExponent = 4
    
    // Proton.
    descriptor.netCharges = [+1]
    descriptor.netSpinPolarizations = [0]
    let proton = Ansatz(descriptor: descriptor)
    XCTAssertEqual(proton.spinDownWaveFunctions.count, 0)
    XCTAssertEqual(proton.spinNeutralWaveFunctions.count, 0)
    XCTAssertEqual(proton.spinUpWaveFunctions.count, 0)
    
    // Neutral hydrogen atom (positive spin).
    descriptor.netCharges = [0]
    descriptor.netSpinPolarizations = [1]
    let hydrogenUp = Ansatz(descriptor: descriptor)
    XCTAssertEqual(hydrogenUp.spinDownWaveFunctions.count, 0)
    XCTAssertEqual(hydrogenUp.spinNeutralWaveFunctions.count, 0)
    XCTAssertEqual(hydrogenUp.spinUpWaveFunctions.count, 1)
    Self.checkFragments(hydrogenUp.spinUpWaveFunctions[0], 1000)
    XCTAssertEqual(1, Self.queryRadius(
      waveFunction: hydrogenUp.spinUpWaveFunctions[0]), accuracy: 0.02)
    
    // Neutral hydrogen atom (negative spin).
    descriptor.netCharges = [0]
    descriptor.netSpinPolarizations = [-1]
    let hydrogenDown = Ansatz(descriptor: descriptor)
    XCTAssertEqual(hydrogenDown.spinDownWaveFunctions.count, 1)
    XCTAssertEqual(hydrogenDown.spinNeutralWaveFunctions.count, 0)
    XCTAssertEqual(hydrogenDown.spinUpWaveFunctions.count, 0)
    Self.checkFragments(hydrogenDown.spinDownWaveFunctions[0], 1000)
    XCTAssertEqual(1, Self.queryRadius(
      waveFunction: hydrogenDown.spinDownWaveFunctions[0]), accuracy: 0.02)
    
    // Hydride anion (singlet state).
    descriptor.netCharges = [-1]
    descriptor.netSpinPolarizations = [0]
    let hydrideSinglet = Ansatz(descriptor: descriptor)
    XCTAssertEqual(hydrideSinglet.spinDownWaveFunctions.count, 0)
    XCTAssertEqual(hydrideSinglet.spinNeutralWaveFunctions.count, 1)
    XCTAssertEqual(hydrideSinglet.spinUpWaveFunctions.count, 0)
    Self.checkFragments(hydrideSinglet.spinNeutralWaveFunctions[0], 1000)
    XCTAssertEqual(1, Self.queryRadius(
      waveFunction: hydrideSinglet.spinNeutralWaveFunctions[0]), accuracy: 0.02)
    
    // Hydride anion (triplet state).
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
    var descriptor = AnsatzDescriptor()
    descriptor.atomicNumbers = [3]
    descriptor.fragmentCount = 1000
    descriptor.positions = [.zero]
    descriptor.sizeExponent = 4
    
    // Metallic ion.
    descriptor.netCharges = [+1]
    descriptor.netSpinPolarizations = [0]
    let lithiumIon = Ansatz(descriptor: descriptor)
    XCTAssertEqual(lithiumIon.spinDownWaveFunctions.count, 0)
    XCTAssertEqual(lithiumIon.spinNeutralWaveFunctions.count, 1)
    XCTAssertEqual(lithiumIon.spinUpWaveFunctions.count, 0)
    Self.checkFragments(lithiumIon.spinNeutralWaveFunctions[0], 1000)
    XCTAssertEqual(0.414, Self.queryRadius(
      waveFunction: lithiumIon.spinNeutralWaveFunctions[0]), accuracy: 0.02)
    
    // Neutral atom.
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
    
    // Spin-3/2 atom.
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
    var descriptor = AnsatzDescriptor()
    descriptor.atomicNumbers = [7, 7]
    descriptor.fragmentCount = 1000
    descriptor.positions = [SIMD3(-1, 0, 0), SIMD3(1, 0, 0)]
    descriptor.sizeExponent = 4
    
    // N+ and N-.
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
    
    // N(3+) anions.
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
    
    // Triplet state.
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
    
    // ROHF simulation that forms 3 covalent bonds.
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
    var descriptor = AnsatzDescriptor()
    descriptor.atomicNumbers = [24]
    descriptor.fragmentCount = 1000
    descriptor.positions = [.zero]
    descriptor.sizeExponent = 4
    
    // 3+ oxidation state, with minimal spin polarization.
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
    
    // Aufbau polarized atom.
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
    
    // Actual polarized atom.
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
  func testGroupIV() throws {
    var descriptor = AnsatzDescriptor()
    descriptor.fragmentCount = 1000
    descriptor.netCharges = [0]
    descriptor.netSpinPolarizations = [2]
    descriptor.positions = [.zero]
    descriptor.sizeExponent = 4
    
    // Carbon diradical atom.
    /*
     Converged (1,000,000 fragments)
     expectation radius: 0.18462892
     expectation radius: 0.99920744
     */
    descriptor.atomicNumbers = [6]
    let carbon = Ansatz(descriptor: descriptor)
    XCTAssertEqual(carbon.spinDownWaveFunctions.count, 0)
    XCTAssertEqual(carbon.spinNeutralWaveFunctions.count, 2)
    XCTAssertEqual(carbon.spinUpWaveFunctions.count, 2)
    for i in 0..<2 {
      Self.checkFragments(carbon.spinNeutralWaveFunctions[i], 1000)
    }
    Self.checkFragments(carbon.spinUpWaveFunctions[0], 1000)
    Self.checkFragments(carbon.spinUpWaveFunctions[1], 1000)
    XCTAssertEqual(1 / 5.414, Self.queryRadius(
      waveFunction: carbon.spinNeutralWaveFunctions[0]), accuracy: 0.01)
    XCTAssertEqual(1, Self.queryRadius(
      waveFunction: carbon.spinNeutralWaveFunctions[1]), accuracy: 0.04)
    
    // Silicon diradical atom.
    /*
     Converged (1,000,000 fragments)
     expectation radius: 0.07451529
     expectation radius: 0.29279003
     expectation radius: 0.24394919
     expectation radius: 0.24394919
     expectation radius: 0.24394919
     expectation radius: 1.4992918
     */
    descriptor.atomicNumbers = [14]
    let silicon = Ansatz(descriptor: descriptor)
    XCTAssertEqual(silicon.spinDownWaveFunctions.count, 0)
    XCTAssertEqual(silicon.spinNeutralWaveFunctions.count, 6)
    XCTAssertEqual(silicon.spinUpWaveFunctions.count, 2)
    for i in 0..<6 {
      Self.checkFragments(silicon.spinNeutralWaveFunctions[i], 1000)
    }
    Self.checkFragments(silicon.spinUpWaveFunctions[0], 1000)
    Self.checkFragments(silicon.spinUpWaveFunctions[1], 1000)
    XCTAssertEqual(1 / 13.414, Self.queryRadius(
      waveFunction: silicon.spinNeutralWaveFunctions[0]), accuracy: 3e-3)
    XCTAssertEqual(2 / 6.828, Self.queryRadius(
      waveFunction: silicon.spinNeutralWaveFunctions[1]), accuracy: 6e-3)
    XCTAssertEqual(2 / 6.828 * 5 / 6, Self.queryRadius(
      waveFunction: silicon.spinNeutralWaveFunctions[2]), accuracy: 6e-3)
    XCTAssertEqual(2 / 6.828 * 5 / 6, Self.queryRadius(
      waveFunction: silicon.spinNeutralWaveFunctions[3]), accuracy: 6e-3)
    XCTAssertEqual(2 / 6.828 * 5 / 6, Self.queryRadius(
      waveFunction: silicon.spinNeutralWaveFunctions[4]), accuracy: 6e-3)
    XCTAssertEqual(1.5, Self.queryRadius(
      waveFunction: silicon.spinNeutralWaveFunctions[5]), accuracy: 0.04)
    
    // Germanium diradical atom.
    /*
     Converged (1,000,000 fragments)
     expectation radius: 0.031816244
     expectation radius: 0.080525205
     expectation radius: 0.06710341
     expectation radius: 0.06710341
     expectation radius: 0.06710341
     expectation radius: 0.36382768
     expectation radius: 0.3367374
     expectation radius: 0.3367374
     expectation radius: 0.3367374
     expectation radius: 0.28296816
     expectation radius: 0.28296816
     expectation radius: 0.2829506
     expectation radius: 0.28296816
     expectation radius: 0.28294447
     expectation radius: 1.9990631
     */
    descriptor.atomicNumbers = [32]
    let germanium = Ansatz(descriptor: descriptor)
    XCTAssertEqual(germanium.spinDownWaveFunctions.count, 0)
    XCTAssertEqual(germanium.spinNeutralWaveFunctions.count, 15)
    XCTAssertEqual(germanium.spinUpWaveFunctions.count, 2)
    for i in 0..<15 {
      Self.checkFragments(germanium.spinNeutralWaveFunctions[i], 1000)
    }
    Self.checkFragments(germanium.spinUpWaveFunctions[0], 1000)
    Self.checkFragments(germanium.spinUpWaveFunctions[1], 1000)
    XCTAssertEqual(1 / 31.414, Self.queryRadius(
      waveFunction: germanium.spinNeutralWaveFunctions[0]), accuracy: 1e-3)
    XCTAssertEqual(2 / 24.828, Self.queryRadius(
      waveFunction: germanium.spinNeutralWaveFunctions[1]), accuracy: 4e-3)
    XCTAssertEqual(2 / 24.828 * 5 / 6, Self.queryRadius(
      waveFunction: germanium.spinNeutralWaveFunctions[2]), accuracy: 4e-3)
    XCTAssertEqual(2 / 24.828 * 5 / 6, Self.queryRadius(
      waveFunction: germanium.spinNeutralWaveFunctions[3]), accuracy: 4e-3)
    XCTAssertEqual(2 / 24.828 * 5 / 6, Self.queryRadius(
      waveFunction: germanium.spinNeutralWaveFunctions[4]), accuracy: 4e-3)
    XCTAssertEqual(0.36382768, Self.queryRadius(
      waveFunction: germanium.spinNeutralWaveFunctions[5]), accuracy: 0.02)
    XCTAssertEqual(0.3367374, Self.queryRadius(
      waveFunction: germanium.spinNeutralWaveFunctions[6]), accuracy: 0.02)
    XCTAssertEqual(0.3367374, Self.queryRadius(
      waveFunction: germanium.spinNeutralWaveFunctions[7]), accuracy: 0.02)
    XCTAssertEqual(0.3367374, Self.queryRadius(
      waveFunction: germanium.spinNeutralWaveFunctions[8]), accuracy: 0.02)
    XCTAssertEqual(0.28296816, Self.queryRadius(
      waveFunction: germanium.spinNeutralWaveFunctions[9]), accuracy: 0.02)
    XCTAssertEqual(0.28296816, Self.queryRadius(
      waveFunction: germanium.spinNeutralWaveFunctions[10]), accuracy: 0.02)
    XCTAssertEqual(0.2829506, Self.queryRadius(
      waveFunction: germanium.spinNeutralWaveFunctions[11]), accuracy: 0.02)
    XCTAssertEqual(0.28296816, Self.queryRadius(
      waveFunction: germanium.spinNeutralWaveFunctions[12]), accuracy: 0.02)
    XCTAssertEqual(0.28294447, Self.queryRadius(
      waveFunction: germanium.spinNeutralWaveFunctions[13]), accuracy: 0.02)
    XCTAssertEqual(2, Self.queryRadius(
      waveFunction: germanium.spinNeutralWaveFunctions[14]), accuracy: 0.04)
    
    // Tin diradical atom.
    /*
     Converged (1,000,000 fragments)
     expectation radius: 0.020224972
     expectation radius: 0.046682257
     expectation radius: 0.03889774
     expectation radius: 0.03889774
     expectation radius: 0.03889774
     expectation radius: 0.114265166
     expectation radius: 0.105806686
     expectation radius: 0.105806686
     expectation radius: 0.105806686
     expectation radius: 0.08886481
     expectation radius: 0.08886481
     expectation radius: 0.08887705
     expectation radius: 0.08886481
     expectation radius: 0.088864505
     expectation radius: 0.48505315
     expectation radius: 0.46485677
     expectation radius: 0.46485677
     expectation radius: 0.46485677
     expectation radius: 0.42438996
     expectation radius: 0.42438996
     expectation radius: 0.42439926
     expectation radius: 0.42438996
     expectation radius: 0.42441788
     expectation radius: 2.4987237
     */
    descriptor.atomicNumbers = [50]
    let tin = Ansatz(descriptor: descriptor)
    XCTAssertEqual(tin.spinDownWaveFunctions.count, 0)
    XCTAssertEqual(tin.spinNeutralWaveFunctions.count, 24)
    XCTAssertEqual(tin.spinUpWaveFunctions.count, 2)
    for i in 0..<24 {
      Self.checkFragments(tin.spinNeutralWaveFunctions[i], 1000)
    }
    Self.checkFragments(tin.spinUpWaveFunctions[0], 1000)
    Self.checkFragments(tin.spinUpWaveFunctions[1], 1000)
    XCTAssertEqual(1 / 49.414, Self.queryRadius(
      waveFunction: tin.spinNeutralWaveFunctions[0]), accuracy: 5e-4)
    XCTAssertEqual(2 / 42.828, Self.queryRadius(
      waveFunction: tin.spinNeutralWaveFunctions[1]), accuracy: 1e-3)
    XCTAssertEqual(2 / 42.828 * 5 / 6, Self.queryRadius(
      waveFunction: tin.spinNeutralWaveFunctions[2]), accuracy: 1e-3)
    XCTAssertEqual(2 / 42.828 * 5 / 6, Self.queryRadius(
      waveFunction: tin.spinNeutralWaveFunctions[3]), accuracy: 1e-3)
    XCTAssertEqual(2 / 42.828 * 5 / 6, Self.queryRadius(
      waveFunction: tin.spinNeutralWaveFunctions[4]), accuracy: 1e-3)
    XCTAssertEqual(0.114265166, Self.queryRadius(
      waveFunction: tin.spinNeutralWaveFunctions[5]), accuracy: 6e-3)
    XCTAssertEqual(0.105806686, Self.queryRadius(
      waveFunction: tin.spinNeutralWaveFunctions[6]), accuracy: 6e-3)
    XCTAssertEqual(0.105806686, Self.queryRadius(
      waveFunction: tin.spinNeutralWaveFunctions[7]), accuracy: 6e-3)
    XCTAssertEqual(0.105806686, Self.queryRadius(
      waveFunction: tin.spinNeutralWaveFunctions[8]), accuracy: 6e-3)
    XCTAssertEqual(0.08886481, Self.queryRadius(
      waveFunction: tin.spinNeutralWaveFunctions[9]), accuracy: 6e-3)
    XCTAssertEqual(0.08886481, Self.queryRadius(
      waveFunction: tin.spinNeutralWaveFunctions[10]), accuracy: 6e-3)
    XCTAssertEqual(0.08887705, Self.queryRadius(
      waveFunction: tin.spinNeutralWaveFunctions[11]), accuracy: 6e-3)
    XCTAssertEqual(0.08886481, Self.queryRadius(
      waveFunction: tin.spinNeutralWaveFunctions[12]), accuracy: 6e-3)
    XCTAssertEqual(0.088864505, Self.queryRadius(
      waveFunction: tin.spinNeutralWaveFunctions[13]), accuracy: 6e-3)
    XCTAssertEqual(2.5, Self.queryRadius(
      waveFunction: tin.spinNeutralWaveFunctions[23]), accuracy: 0.06)
    
    // Lead diradical atom.
    /*
     Converged (1,000,000 fragments)
     expectation radius: 0.01227764
     expectation radius: 0.026715986
     expectation radius: 0.022263272
     expectation radius: 0.022263272
     expectation radius: 0.022263272
     expectation radius: 0.051483016
     expectation radius: 0.04767065
     expectation radius: 0.04767065
     expectation radius: 0.04767065
     expectation radius: 0.0400402
     expectation radius: 0.0400402
     expectation radius: 0.040046137
     expectation radius: 0.0400402
     expectation radius: 0.040042643
     expectation radius: 0.14456497
     expectation radius: 0.13852988
     expectation radius: 0.13852988
     expectation radius: 0.13852988
     expectation radius: 0.1264886
     expectation radius: 0.1264886
     expectation radius: 0.12646745
     expectation radius: 0.1264886
     expectation radius: 0.12648587
     expectation radius: 0.1084131
     expectation radius: 0.10839141
     expectation radius: 0.10841278
     expectation radius: 0.10840806
     expectation radius: 0.10841278
     expectation radius: 0.108416915
     expectation radius: 0.1084131
     expectation radius: 0.60633916
     expectation radius: 0.5901139
     expectation radius: 0.5901139
     expectation radius: 0.5901139
     expectation radius: 0.5577721
     expectation radius: 0.5577721
     expectation radius: 0.55763835
     expectation radius: 0.5577721
     expectation radius: 0.55775756
     expectation radius: 2.9965723
     */
    descriptor.atomicNumbers = [82]
    let lead = Ansatz(descriptor: descriptor)
    XCTAssertEqual(lead.spinDownWaveFunctions.count, 0)
    XCTAssertEqual(lead.spinNeutralWaveFunctions.count, 40)
    XCTAssertEqual(lead.spinUpWaveFunctions.count, 2)
    for i in 0..<40 {
      Self.checkFragments(lead.spinNeutralWaveFunctions[i], 1000, 5000)
    }
    Self.checkFragments(lead.spinUpWaveFunctions[0], 1000)
    Self.checkFragments(lead.spinUpWaveFunctions[1], 1000)
    XCTAssertEqual(1 / 81.414, Self.queryRadius(
      waveFunction: lead.spinNeutralWaveFunctions[0]), accuracy: 3e-4)
    XCTAssertEqual(2 / 74.828, Self.queryRadius(
      waveFunction: lead.spinNeutralWaveFunctions[1]), accuracy: 1e-3)
    XCTAssertEqual(2 / 74.828 * 5 / 6, Self.queryRadius(
      waveFunction: lead.spinNeutralWaveFunctions[2]), accuracy: 1e-3)
    XCTAssertEqual(2 / 74.828 * 5 / 6, Self.queryRadius(
      waveFunction: lead.spinNeutralWaveFunctions[3]), accuracy: 1e-3)
    XCTAssertEqual(2 / 74.828 * 5 / 6, Self.queryRadius(
      waveFunction: lead.spinNeutralWaveFunctions[4]), accuracy: 1e-3)
    XCTAssertEqual(0.051483016, Self.queryRadius(
      waveFunction: lead.spinNeutralWaveFunctions[5]), accuracy: 3e-3)
    XCTAssertEqual(0.04767065, Self.queryRadius(
      waveFunction: lead.spinNeutralWaveFunctions[6]), accuracy: 3e-3)
    XCTAssertEqual(0.04767065, Self.queryRadius(
      waveFunction: lead.spinNeutralWaveFunctions[7]), accuracy: 3e-3)
    XCTAssertEqual(0.04767065, Self.queryRadius(
      waveFunction: lead.spinNeutralWaveFunctions[8]), accuracy: 3e-3)
    XCTAssertEqual(0.0400402, Self.queryRadius(
      waveFunction: lead.spinNeutralWaveFunctions[9]), accuracy: 3e-3)
    XCTAssertEqual(0.0400402, Self.queryRadius(
      waveFunction: lead.spinNeutralWaveFunctions[10]), accuracy: 3e-3)
    XCTAssertEqual(0.040046137, Self.queryRadius(
      waveFunction: lead.spinNeutralWaveFunctions[11]), accuracy: 3e-3)
    XCTAssertEqual(0.0400402, Self.queryRadius(
      waveFunction: lead.spinNeutralWaveFunctions[12]), accuracy: 3e-3)
    XCTAssertEqual(0.040042643, Self.queryRadius(
      waveFunction: lead.spinNeutralWaveFunctions[13]), accuracy: 3e-3)
    XCTAssertEqual(3, Self.queryRadius(
      waveFunction: lead.spinNeutralWaveFunctions[39]), accuracy: 0.06)
    
    // Flerovium ion.
    //
    // Challenging test case for the performance of this framework: process
    // the superheavy element flerovium. This triggers the most severe edge
    // cases with wavefunctions fluctuating over large length scales.
    //
    // We need 5000 fragments to have good expectation radii for flerovium.
    /*
     Converged (10,000 fragments)
     expectation radius: 0.008773583
     expectation radius: 0.018610345
     expectation radius: 0.0154791875
     expectation radius: 0.0154791875
     expectation radius: 0.0154791875
     expectation radius: 0.03310816
     expectation radius: 0.03059156
     expectation radius: 0.03059156
     expectation radius: 0.03059156
     expectation radius: 0.02567632
     expectation radius: 0.02567632
     expectation radius: 0.025730968
     expectation radius: 0.02567632
     expectation radius: 0.025719682
     expectation radius: 0.049063906
     expectation radius: 0.04700585
     expectation radius: 0.04700585
     expectation radius: 0.04700585
     
     Converged (100,000 fragments)
     expectation radius: 0.008803415
     expectation radius: 0.018659757
     expectation radius: 0.015560842
     expectation radius: 0.015560842
     expectation radius: 0.015560842
     expectation radius: 0.033196677
     expectation radius: 0.03072553
     expectation radius: 0.03072553
     expectation radius: 0.03072553
     expectation radius: 0.025792854
     expectation radius: 0.025792854
     expectation radius: 0.02580407
     expectation radius: 0.025792854
     expectation radius: 0.025808508
     expectation radius: 0.049400333
     expectation radius: 0.047327485
     expectation radius: 0.047327485
     expectation radius: 0.047327485
     
     Converged (1,000,000 fragments)
     expectation radius: 0.008813468
     expectation radius: 0.018714862
     expectation radius: 0.015593298
     expectation radius: 0.015593298
     expectation radius: 0.015593298
     expectation radius: 0.033228643
     expectation radius: 0.030767119
     expectation radius: 0.030767119
     expectation radius: 0.030767119
     expectation radius: 0.02584201
     expectation radius: 0.02584201
     expectation radius: 0.025841631
     expectation radius: 0.02584201
     expectation radius: 0.02584289
     expectation radius: 0.049465608
     expectation radius: 0.047402453
     expectation radius: 0.047402453
     expectation radius: 0.047402453
     */
    descriptor.fragmentCount = 5000
    descriptor.atomicNumbers = [114]
    descriptor.netCharges = [78]
    descriptor.netSpinPolarizations = [0]
    let flerovium = Ansatz(descriptor: descriptor)
    XCTAssertEqual(flerovium.spinDownWaveFunctions.count, 0)
    XCTAssertEqual(flerovium.spinNeutralWaveFunctions.count, 18)
    XCTAssertEqual(flerovium.spinUpWaveFunctions.count, 0)
    for i in 0..<18 {
      Self.checkFragments(flerovium.spinNeutralWaveFunctions[i], 5000, 30000)
    }
    XCTAssertEqual(1 / 113.414, Self.queryRadius(
      waveFunction: flerovium.spinNeutralWaveFunctions[0]), accuracy: 1e-4)
    XCTAssertEqual(2 / 106.828, Self.queryRadius(
      waveFunction: flerovium.spinNeutralWaveFunctions[1]), accuracy: 3e-4)
    XCTAssertEqual(2 / 106.828 * 5 / 6, Self.queryRadius(
      waveFunction: flerovium.spinNeutralWaveFunctions[2]), accuracy: 3e-4)
    XCTAssertEqual(2 / 106.828 * 5 / 6, Self.queryRadius(
      waveFunction: flerovium.spinNeutralWaveFunctions[3]), accuracy: 3e-4)
    XCTAssertEqual(2 / 106.828 * 5 / 6, Self.queryRadius(
      waveFunction: flerovium.spinNeutralWaveFunctions[4]), accuracy: 3e-4)
    XCTAssertEqual(0.033228643, Self.queryRadius(
      waveFunction: flerovium.spinNeutralWaveFunctions[5]), accuracy: 1e-3)
    XCTAssertEqual(0.030767119, Self.queryRadius(
      waveFunction: flerovium.spinNeutralWaveFunctions[6]), accuracy: 1e-3)
    XCTAssertEqual(0.030767119, Self.queryRadius(
      waveFunction: flerovium.spinNeutralWaveFunctions[7]), accuracy: 1e-3)
    XCTAssertEqual(0.030767119, Self.queryRadius(
      waveFunction: flerovium.spinNeutralWaveFunctions[8]), accuracy: 1e-3)
    XCTAssertEqual(0.02584201, Self.queryRadius(
      waveFunction: flerovium.spinNeutralWaveFunctions[9]), accuracy: 1e-3)
    XCTAssertEqual(0.02584201, Self.queryRadius(
      waveFunction: flerovium.spinNeutralWaveFunctions[10]), accuracy: 1e-3)
    XCTAssertEqual(0.025841631, Self.queryRadius(
      waveFunction: flerovium.spinNeutralWaveFunctions[11]), accuracy: 1e-3)
    XCTAssertEqual(0.02584201, Self.queryRadius(
      waveFunction: flerovium.spinNeutralWaveFunctions[12]), accuracy: 1e-3)
    XCTAssertEqual(0.02584289, Self.queryRadius(
      waveFunction: flerovium.spinNeutralWaveFunctions[13]), accuracy: 1e-3)
    XCTAssertEqual(0.049465608, Self.queryRadius(
      waveFunction: flerovium.spinNeutralWaveFunctions[14]), accuracy: 1e-3)
    XCTAssertEqual(0.047402453, Self.queryRadius(
      waveFunction: flerovium.spinNeutralWaveFunctions[15]), accuracy: 1e-3)
    XCTAssertEqual(0.047402453, Self.queryRadius(
      waveFunction: flerovium.spinNeutralWaveFunctions[16]), accuracy: 1e-3)
    XCTAssertEqual(0.047402453, Self.queryRadius(
      waveFunction: flerovium.spinNeutralWaveFunctions[17]), accuracy: 1e-3)
  }
  
  // Next, optimize the code for generating the initial guess.
}
