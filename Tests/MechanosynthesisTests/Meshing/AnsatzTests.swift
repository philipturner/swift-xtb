import XCTest
import Mechanosynthesis
import Numerics

final class AnsatzTests: XCTestCase {
  static func checkFragments(
    _ orbital: HydrogenicOrbital,
    _ expectedCount: Int
  ) {
    var octreeFragmentCount = 0
    for node in orbital.octree.nodes {
      let mask8 = node.branchesMask & SIMD8(repeating: 128)
      let mask64 = unsafeBitCast(mask8, to: UInt64.self)
      octreeFragmentCount += 8 * mask64.nonzeroBitCount
    }
    XCTAssertGreaterThanOrEqual(octreeFragmentCount, expectedCount)
    XCTAssertLessThan(octreeFragmentCount, expectedCount * 2)
  }
  
  static func queryRadius(
    orbital: HydrogenicOrbital,
    nucleusPosition: SIMD3<Float> = .zero
  ) -> Float {
    var sum: Double = .zero
    for nodeID in orbital.octree.nodes.indices {
      let node = orbital.octree.nodes[nodeID]
      
      // Lookup table for child nodes.
      var lx = SIMD8<Float>(0, 1, 0, 1, 0, 1, 0, 1) * 0.5 - 0.25
      var ly = SIMD8<Float>(0, 0, 1, 1, 0, 0, 1, 1) * 0.5 - 0.25
      var lz = SIMD8<Float>(0, 0, 0, 0, 1, 1, 1, 1) * 0.5 - 0.25
      let centerDelta = node.center - nucleusPosition
      lx = lx * node.spacing + centerDelta.x
      ly = ly * node.spacing + centerDelta.y
      lz = lz * node.spacing + centerDelta.z
      
      // Evaluate the radius at 2x the resolution of octree nodes.
      for branchID in 0..<8
      where node.branchesMask[branchID] == UInt8.max {
        var x = SIMD8<Float>(0, 1, 0, 1, 0, 1, 0, 1) * 0.5 - 0.25
        var y = SIMD8<Float>(0, 0, 1, 1, 0, 0, 1, 1) * 0.5 - 0.25
        var z = SIMD8<Float>(0, 0, 0, 0, 1, 1, 1, 1) * 0.5 - 0.25
        x = x * node.spacing / 2 + lx[branchID]
        y = y * node.spacing / 2 + ly[branchID]
        z = z * node.spacing / 2 + lz[branchID]
        
        let Ψ = orbital.basisFunction.amplitude(x: x, y: y, z: z)
        let d3r = node.spacing * node.spacing * node.spacing / 64
        let r = (x * x + y * y + z * z).squareRoot()
        let ΨrΨ = Ψ * r * Ψ * d3r
        sum += Double(ΨrΨ.sum())
      }
    }
    
    return Float(sum) * 2 / 3
  }
  
  func testHydrogen() throws {
    var descriptor = AnsatzDescriptor()
    descriptor.atomicNumber = 1
    descriptor.fragmentCount = 5000
    descriptor.position = .zero
    descriptor.sizeExponent = 4
    
    // Proton.
    descriptor.netCharge = +1
    descriptor.netSpin = 0
    let proton = Ansatz(descriptor: descriptor)
    XCTAssertEqual(proton.spinDownOrbitals.count, 0)
    XCTAssertEqual(proton.spinNeutralOrbitals.count, 0)
    XCTAssertEqual(proton.spinUpOrbitals.count, 0)
    
    // Neutral hydrogen atom (positive spin).
    descriptor.netCharge = 0
    descriptor.netSpin = 0.5
    let hydrogenUp = Ansatz(descriptor: descriptor)
    XCTAssertEqual(hydrogenUp.spinDownOrbitals.count, 0)
    XCTAssertEqual(hydrogenUp.spinNeutralOrbitals.count, 0)
    XCTAssertEqual(hydrogenUp.spinUpOrbitals.count, 1)
    Self.checkFragments(hydrogenUp.spinUpOrbitals[0], 5000)
    XCTAssertEqual(1, Self.queryRadius(
      orbital: hydrogenUp.spinUpOrbitals[0]), accuracy: 0.02)
    
    // Neutral hydrogen atom (negative spin).
    descriptor.netCharge = 0
    descriptor.netSpin = -0.5
    let hydrogenDown = Ansatz(descriptor: descriptor)
    XCTAssertEqual(hydrogenDown.spinDownOrbitals.count, 1)
    XCTAssertEqual(hydrogenDown.spinNeutralOrbitals.count, 0)
    XCTAssertEqual(hydrogenDown.spinUpOrbitals.count, 0)
    Self.checkFragments(hydrogenDown.spinDownOrbitals[0], 5000)
    XCTAssertEqual(1, Self.queryRadius(
      orbital: hydrogenDown.spinDownOrbitals[0]), accuracy: 0.02)
    
    // Hydride anion (singlet state).
    descriptor.netCharge = -1
    descriptor.netSpin = 0
    let hydrideSinglet = Ansatz(descriptor: descriptor)
    XCTAssertEqual(hydrideSinglet.spinDownOrbitals.count, 0)
    XCTAssertEqual(hydrideSinglet.spinNeutralOrbitals.count, 1)
    XCTAssertEqual(hydrideSinglet.spinUpOrbitals.count, 0)
    Self.checkFragments(hydrideSinglet.spinNeutralOrbitals[0], 5000)
    XCTAssertEqual(1, Self.queryRadius(
      orbital: hydrideSinglet.spinNeutralOrbitals[0]), accuracy: 0.02)
    
    // Hydride anion (triplet state).
    descriptor.atomicNumber = 1
    descriptor.netCharge = -1
    descriptor.netSpin = 1
    let hydrideTriplet = Ansatz(descriptor: descriptor)
    XCTAssertEqual(hydrideTriplet.spinDownOrbitals.count, 0)
    XCTAssertEqual(hydrideTriplet.spinNeutralOrbitals.count, 0)
    XCTAssertEqual(hydrideTriplet.spinUpOrbitals.count, 2)
    Self.checkFragments(hydrideTriplet.spinUpOrbitals[0], 5000)
    Self.checkFragments(hydrideTriplet.spinUpOrbitals[1], 5000)
    XCTAssertEqual(1, Self.queryRadius(
      orbital: hydrideTriplet.spinUpOrbitals[0]), accuracy: 0.02)
    XCTAssertEqual(2, Self.queryRadius(
      orbital: hydrideTriplet.spinUpOrbitals[1]), accuracy: 0.08)
  }
  
  func testLithium() throws {
    var descriptor = AnsatzDescriptor()
    descriptor.atomicNumber = 3
    descriptor.fragmentCount = 5000
    descriptor.position = .zero
    descriptor.sizeExponent = 4
    
    // Metallic ion.
    descriptor.netCharge = +1
    descriptor.netSpin = 0
    let lithiumIon = Ansatz(descriptor: descriptor)
    XCTAssertEqual(lithiumIon.spinDownOrbitals.count, 0)
    XCTAssertEqual(lithiumIon.spinNeutralOrbitals.count, 1)
    XCTAssertEqual(lithiumIon.spinUpOrbitals.count, 0)
    Self.checkFragments(lithiumIon.spinNeutralOrbitals[0], 5000)
    XCTAssertEqual(0.414, Self.queryRadius(
      orbital: lithiumIon.spinNeutralOrbitals[0]), accuracy: 0.02)
    
    // Neutral atom.
    descriptor.netCharge = 0
    descriptor.netSpin = -0.5
    let lithiumNeutral = Ansatz(descriptor: descriptor)
    XCTAssertEqual(lithiumNeutral.spinDownOrbitals.count, 1)
    XCTAssertEqual(lithiumNeutral.spinNeutralOrbitals.count, 1)
    XCTAssertEqual(lithiumNeutral.spinUpOrbitals.count, 0)
    Self.checkFragments(lithiumNeutral.spinDownOrbitals[0], 5000)
    Self.checkFragments(lithiumNeutral.spinNeutralOrbitals[0], 5000)
    XCTAssertEqual(0.414, Self.queryRadius(
      orbital: lithiumNeutral.spinNeutralOrbitals[0]), accuracy: 0.02)
    XCTAssertEqual(2, Self.queryRadius(
      orbital: lithiumNeutral.spinDownOrbitals[0]), accuracy: 0.08)
    
    // Spin-3/2 atom.
    descriptor.netCharge = 0
    descriptor.netSpin = 1.5
    let lithiumPolarized = Ansatz(descriptor: descriptor)
    XCTAssertEqual(lithiumPolarized.spinDownOrbitals.count, 0)
    XCTAssertEqual(lithiumPolarized.spinNeutralOrbitals.count, 0)
    XCTAssertEqual(lithiumPolarized.spinUpOrbitals.count, 3)
    Self.checkFragments(lithiumPolarized.spinUpOrbitals[0], 5000)
    Self.checkFragments(lithiumPolarized.spinUpOrbitals[1], 5000)
    Self.checkFragments(lithiumPolarized.spinUpOrbitals[2], 5000)
    XCTAssertEqual(0.333, Self.queryRadius(
      orbital: lithiumPolarized.spinUpOrbitals[0]), accuracy: 0.02)
    XCTAssertEqual(1.414, Self.queryRadius(
      orbital: lithiumPolarized.spinUpOrbitals[1]), accuracy: 0.03)
    XCTAssertEqual(1.414 * 5 / 6, Self.queryRadius(
      orbital: lithiumPolarized.spinUpOrbitals[2]), accuracy: 0.03)
  }
  
  // Test chromium.
  // - Test the limitation of GFN-xTB, which doesn't support the true spin
  //   polarization of chromium d-orbitals. Instead, it does 3d^4, 4s^2.
  func testChromium() throws {
    var descriptor = AnsatzDescriptor()
    descriptor.atomicNumber = 24
    descriptor.fragmentCount = 5000
    descriptor.position = .zero
    descriptor.sizeExponent = 4
    
    // 3+ oxidation state, with minimal spin polarization.
    descriptor.netCharge = +3
    descriptor.netSpin = 0.5
    let chromiumIon = Ansatz(descriptor: descriptor)
    XCTAssertEqual(chromiumIon.spinDownOrbitals.count, 0)
    XCTAssertEqual(chromiumIon.spinNeutralOrbitals.count, 10)
    XCTAssertEqual(chromiumIon.spinUpOrbitals.count, 1)
    for i in 0..<10 {
      Self.checkFragments(chromiumIon.spinNeutralOrbitals[i], 5000)
    }
    Self.checkFragments(chromiumIon.spinUpOrbitals[0], 5000)
    XCTAssertEqual(1 / 23.414, Self.queryRadius(
      orbital: chromiumIon.spinNeutralOrbitals[0]), accuracy: 2e-3)
    XCTAssertEqual(2 / 16.828, Self.queryRadius(
      orbital: chromiumIon.spinNeutralOrbitals[1]), accuracy: 4e-3)
    XCTAssertEqual(2 / 16.828 * 5 / 6, Self.queryRadius(
      orbital: chromiumIon.spinNeutralOrbitals[2]), accuracy: 4e-3)
    XCTAssertEqual(2 / 16.828 * 5 / 6, Self.queryRadius(
      orbital: chromiumIon.spinNeutralOrbitals[3]), accuracy: 4e-3)
    XCTAssertEqual(2 / 16.828 * 5 / 6, Self.queryRadius(
      orbital: chromiumIon.spinNeutralOrbitals[4]), accuracy: 4e-3)
    XCTAssertEqual(0.375, Self.queryRadius(
      orbital: chromiumIon.spinNeutralOrbitals[5]), accuracy: 0.01)
    XCTAssertEqual(0.375 * 12 / 13, Self.queryRadius(
      orbital: chromiumIon.spinNeutralOrbitals[6]), accuracy: 0.01)
    XCTAssertEqual(0.375 * 12 / 13, Self.queryRadius(
      orbital: chromiumIon.spinNeutralOrbitals[7]), accuracy: 0.01)
    XCTAssertEqual(0.375 * 12 / 13, Self.queryRadius(
      orbital: chromiumIon.spinNeutralOrbitals[8]), accuracy: 0.01)
    XCTAssertEqual(0.2906231, Self.queryRadius(
      orbital: chromiumIon.spinUpOrbitals[0]), accuracy: 0.01)
    XCTAssertEqual(4 / 4.414, Self.queryRadius(
      orbital: chromiumIon.spinNeutralOrbitals[9]), accuracy: 0.04)
    
    // Aufbau polarized atom.
    descriptor.netCharge = 0
    descriptor.netSpin = 2
    let chromiumAufbau = Ansatz(descriptor: descriptor)
    XCTAssertEqual(chromiumAufbau.spinDownOrbitals.count, 0)
    XCTAssertEqual(chromiumAufbau.spinNeutralOrbitals.count, 10)
    XCTAssertEqual(chromiumAufbau.spinUpOrbitals.count, 4)
    for i in 0..<10 {
      Self.checkFragments(chromiumAufbau.spinNeutralOrbitals[i], 5000)
    }
    for i in 0..<4 {
      Self.checkFragments(chromiumAufbau.spinUpOrbitals[i], 5000)
    }
    
    // Actual polarized atom.
    descriptor.netCharge = 0
    descriptor.netSpin = 3
    let chromiumActual = Ansatz(descriptor: descriptor)
    XCTAssertEqual(chromiumActual.spinDownOrbitals.count, 0)
    XCTAssertEqual(chromiumActual.spinNeutralOrbitals.count, 9)
    XCTAssertEqual(chromiumActual.spinUpOrbitals.count, 6)
    for i in 0..<9 {
      Self.checkFragments(chromiumActual.spinNeutralOrbitals[i], 5000)
    }
    for i in 0..<6 {
      Self.checkFragments(chromiumActual.spinUpOrbitals[i], 5000)
    }
  }
  
  // Test group (IV) atoms.
  func testGroupIV() throws {
    var descriptor = AnsatzDescriptor()
    descriptor.fragmentCount = 5000
    descriptor.netCharge = 0
    descriptor.netSpin = 1
    descriptor.position = .zero
    descriptor.sizeExponent = 4
    
    // Carbon diradical atom.
    /*
     Converged (1,000,000 fragments)
     expectation radius: 0.18462892
     expectation radius: 0.99920744
     */
    descriptor.atomicNumber = 6
    let carbon = Ansatz(descriptor: descriptor)
    XCTAssertEqual(carbon.spinDownOrbitals.count, 0)
    XCTAssertEqual(carbon.spinNeutralOrbitals.count, 2)
    XCTAssertEqual(carbon.spinUpOrbitals.count, 2)
    for i in 0..<2 {
      Self.checkFragments(carbon.spinNeutralOrbitals[i], 5000)
    }
    Self.checkFragments(carbon.spinUpOrbitals[0], 5000)
    Self.checkFragments(carbon.spinUpOrbitals[1], 5000)
    XCTAssertEqual(1 / 5.414, Self.queryRadius(
      orbital: carbon.spinNeutralOrbitals[0]), accuracy: 0.01)
    XCTAssertEqual(1, Self.queryRadius(
      orbital: carbon.spinNeutralOrbitals[1]), accuracy: 0.04)
    
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
    descriptor.atomicNumber = 14
    let silicon = Ansatz(descriptor: descriptor)
    XCTAssertEqual(silicon.spinDownOrbitals.count, 0)
    XCTAssertEqual(silicon.spinNeutralOrbitals.count, 6)
    XCTAssertEqual(silicon.spinUpOrbitals.count, 2)
    for i in 0..<6 {
      Self.checkFragments(silicon.spinNeutralOrbitals[i], 5000)
    }
    Self.checkFragments(silicon.spinUpOrbitals[0], 5000)
    Self.checkFragments(silicon.spinUpOrbitals[1], 5000)
    XCTAssertEqual(1 / 13.414, Self.queryRadius(
      orbital: silicon.spinNeutralOrbitals[0]), accuracy: 3e-3)
    XCTAssertEqual(2 / 6.828, Self.queryRadius(
      orbital: silicon.spinNeutralOrbitals[1]), accuracy: 6e-3)
    XCTAssertEqual(2 / 6.828 * 5 / 6, Self.queryRadius(
      orbital: silicon.spinNeutralOrbitals[2]), accuracy: 6e-3)
    XCTAssertEqual(2 / 6.828 * 5 / 6, Self.queryRadius(
      orbital: silicon.spinNeutralOrbitals[3]), accuracy: 6e-3)
    XCTAssertEqual(2 / 6.828 * 5 / 6, Self.queryRadius(
      orbital: silicon.spinNeutralOrbitals[4]), accuracy: 6e-3)
    XCTAssertEqual(1.5, Self.queryRadius(
      orbital: silicon.spinNeutralOrbitals[5]), accuracy: 0.04)
    
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
    descriptor.atomicNumber = 32
    let germanium = Ansatz(descriptor: descriptor)
    XCTAssertEqual(germanium.spinDownOrbitals.count, 0)
    XCTAssertEqual(germanium.spinNeutralOrbitals.count, 15)
    XCTAssertEqual(germanium.spinUpOrbitals.count, 2)
    for i in 0..<15 {
      Self.checkFragments(germanium.spinNeutralOrbitals[i], 5000)
    }
    Self.checkFragments(germanium.spinUpOrbitals[0], 5000)
    Self.checkFragments(germanium.spinUpOrbitals[1], 5000)
    XCTAssertEqual(1 / 31.414, Self.queryRadius(
      orbital: germanium.spinNeutralOrbitals[0]), accuracy: 1e-3)
    XCTAssertEqual(2 / 24.828, Self.queryRadius(
      orbital: germanium.spinNeutralOrbitals[1]), accuracy: 4e-3)
    XCTAssertEqual(2 / 24.828 * 5 / 6, Self.queryRadius(
      orbital: germanium.spinNeutralOrbitals[2]), accuracy: 4e-3)
    XCTAssertEqual(2 / 24.828 * 5 / 6, Self.queryRadius(
      orbital: germanium.spinNeutralOrbitals[3]), accuracy: 4e-3)
    XCTAssertEqual(2 / 24.828 * 5 / 6, Self.queryRadius(
      orbital: germanium.spinNeutralOrbitals[4]), accuracy: 4e-3)
    XCTAssertEqual(0.36382768, Self.queryRadius(
      orbital: germanium.spinNeutralOrbitals[5]), accuracy: 0.02)
    XCTAssertEqual(0.3367374, Self.queryRadius(
      orbital: germanium.spinNeutralOrbitals[6]), accuracy: 0.02)
    XCTAssertEqual(0.3367374, Self.queryRadius(
      orbital: germanium.spinNeutralOrbitals[7]), accuracy: 0.02)
    XCTAssertEqual(0.3367374, Self.queryRadius(
      orbital: germanium.spinNeutralOrbitals[8]), accuracy: 0.02)
    XCTAssertEqual(0.28296816, Self.queryRadius(
      orbital: germanium.spinNeutralOrbitals[9]), accuracy: 0.02)
    XCTAssertEqual(0.28296816, Self.queryRadius(
      orbital: germanium.spinNeutralOrbitals[10]), accuracy: 0.02)
    XCTAssertEqual(0.2829506, Self.queryRadius(
      orbital: germanium.spinNeutralOrbitals[11]), accuracy: 0.02)
    XCTAssertEqual(0.28296816, Self.queryRadius(
      orbital: germanium.spinNeutralOrbitals[12]), accuracy: 0.02)
    XCTAssertEqual(0.28294447, Self.queryRadius(
      orbital: germanium.spinNeutralOrbitals[13]), accuracy: 0.02)
    XCTAssertEqual(2, Self.queryRadius(
      orbital: germanium.spinNeutralOrbitals[14]), accuracy: 0.08)
    
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
    descriptor.fragmentCount = 6000
    descriptor.atomicNumber = 50
    let tin = Ansatz(descriptor: descriptor)
    XCTAssertEqual(tin.spinDownOrbitals.count, 0)
    XCTAssertEqual(tin.spinNeutralOrbitals.count, 24)
    XCTAssertEqual(tin.spinUpOrbitals.count, 2)
    for i in 0..<24 {
      Self.checkFragments(tin.spinNeutralOrbitals[i], 6000)
    }
    Self.checkFragments(tin.spinUpOrbitals[0], 6000)
    Self.checkFragments(tin.spinUpOrbitals[1], 6000)
    XCTAssertEqual(1 / 49.414, Self.queryRadius(
      orbital: tin.spinNeutralOrbitals[0]), accuracy: 5e-4)
    XCTAssertEqual(2 / 42.828, Self.queryRadius(
      orbital: tin.spinNeutralOrbitals[1]), accuracy: 1e-3)
    XCTAssertEqual(2 / 42.828 * 5 / 6, Self.queryRadius(
      orbital: tin.spinNeutralOrbitals[2]), accuracy: 1e-3)
    XCTAssertEqual(2 / 42.828 * 5 / 6, Self.queryRadius(
      orbital: tin.spinNeutralOrbitals[3]), accuracy: 1e-3)
    XCTAssertEqual(2 / 42.828 * 5 / 6, Self.queryRadius(
      orbital: tin.spinNeutralOrbitals[4]), accuracy: 1e-3)
    XCTAssertEqual(0.114265166, Self.queryRadius(
      orbital: tin.spinNeutralOrbitals[5]), accuracy: 6e-3)
    XCTAssertEqual(0.105806686, Self.queryRadius(
      orbital: tin.spinNeutralOrbitals[6]), accuracy: 6e-3)
    XCTAssertEqual(0.105806686, Self.queryRadius(
      orbital: tin.spinNeutralOrbitals[7]), accuracy: 6e-3)
    XCTAssertEqual(0.105806686, Self.queryRadius(
      orbital: tin.spinNeutralOrbitals[8]), accuracy: 6e-3)
    XCTAssertEqual(0.08886481, Self.queryRadius(
      orbital: tin.spinNeutralOrbitals[9]), accuracy: 6e-3)
    XCTAssertEqual(0.08886481, Self.queryRadius(
      orbital: tin.spinNeutralOrbitals[10]), accuracy: 6e-3)
    XCTAssertEqual(0.08887705, Self.queryRadius(
      orbital: tin.spinNeutralOrbitals[11]), accuracy: 6e-3)
    XCTAssertEqual(0.08886481, Self.queryRadius(
      orbital: tin.spinNeutralOrbitals[12]), accuracy: 6e-3)
    XCTAssertEqual(0.088864505, Self.queryRadius(
      orbital: tin.spinNeutralOrbitals[13]), accuracy: 6e-3)
    XCTAssertEqual(2.5, Self.queryRadius(
      orbital: tin.spinNeutralOrbitals[23]), accuracy: 0.06)
    
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
    descriptor.fragmentCount = 5000
    descriptor.atomicNumber = 82
    let lead = Ansatz(descriptor: descriptor)
    XCTAssertEqual(lead.spinDownOrbitals.count, 0)
    XCTAssertEqual(lead.spinNeutralOrbitals.count, 40)
    XCTAssertEqual(lead.spinUpOrbitals.count, 2)
    for i in 0..<40 {
      Self.checkFragments(lead.spinNeutralOrbitals[i], 5000)
    }
    Self.checkFragments(lead.spinUpOrbitals[0], 5000)
    Self.checkFragments(lead.spinUpOrbitals[1], 5000)
    XCTAssertEqual(1 / 81.414, Self.queryRadius(
      orbital: lead.spinNeutralOrbitals[0]), accuracy: 3e-4)
    XCTAssertEqual(2 / 74.828, Self.queryRadius(
      orbital: lead.spinNeutralOrbitals[1]), accuracy: 1e-3)
    XCTAssertEqual(2 / 74.828 * 5 / 6, Self.queryRadius(
      orbital: lead.spinNeutralOrbitals[2]), accuracy: 1e-3)
    XCTAssertEqual(2 / 74.828 * 5 / 6, Self.queryRadius(
      orbital: lead.spinNeutralOrbitals[3]), accuracy: 1e-3)
    XCTAssertEqual(2 / 74.828 * 5 / 6, Self.queryRadius(
      orbital: lead.spinNeutralOrbitals[4]), accuracy: 1e-3)
    XCTAssertEqual(0.051483016, Self.queryRadius(
      orbital: lead.spinNeutralOrbitals[5]), accuracy: 3e-3)
    XCTAssertEqual(0.04767065, Self.queryRadius(
      orbital: lead.spinNeutralOrbitals[6]), accuracy: 3e-3)
    XCTAssertEqual(0.04767065, Self.queryRadius(
      orbital: lead.spinNeutralOrbitals[7]), accuracy: 3e-3)
    XCTAssertEqual(0.04767065, Self.queryRadius(
      orbital: lead.spinNeutralOrbitals[8]), accuracy: 3e-3)
    XCTAssertEqual(0.0400402, Self.queryRadius(
      orbital: lead.spinNeutralOrbitals[9]), accuracy: 3e-3)
    XCTAssertEqual(0.0400402, Self.queryRadius(
      orbital: lead.spinNeutralOrbitals[10]), accuracy: 3e-3)
    XCTAssertEqual(0.040046137, Self.queryRadius(
      orbital: lead.spinNeutralOrbitals[11]), accuracy: 3e-3)
    XCTAssertEqual(0.0400402, Self.queryRadius(
      orbital: lead.spinNeutralOrbitals[12]), accuracy: 3e-3)
    XCTAssertEqual(0.040042643, Self.queryRadius(
      orbital: lead.spinNeutralOrbitals[13]), accuracy: 3e-3)
    XCTAssertEqual(3, Self.queryRadius(
      orbital: lead.spinNeutralOrbitals[39]), accuracy: 0.06)
    
    // Flerovium ion.
    //
    // Challenging test case for the performance of this framework: process
    // the superheavy element flerovium. This triggers the most severe edge
    // cases with wavefunctions fluctuating over large length scales.
    //
    // We need 5,000 fragments to have good expectation radii for flerovium.
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
    descriptor.fragmentCount = 12000
    descriptor.atomicNumber = 114
    descriptor.netCharge = 78
    descriptor.netSpin = 0
    let flerovium = Ansatz(descriptor: descriptor)
    XCTAssertEqual(flerovium.spinDownOrbitals.count, 0)
    XCTAssertEqual(flerovium.spinNeutralOrbitals.count, 18)
    XCTAssertEqual(flerovium.spinUpOrbitals.count, 0)
    for i in 0..<18 {
      Self.checkFragments(flerovium.spinNeutralOrbitals[i], 12000)
    }
    XCTAssertEqual(1 / 113.414, Self.queryRadius(
      orbital: flerovium.spinNeutralOrbitals[0]), accuracy: 2e-4)
    XCTAssertEqual(2 / 106.828, Self.queryRadius(
      orbital: flerovium.spinNeutralOrbitals[1]), accuracy: 3e-4)
    XCTAssertEqual(2 / 106.828 * 5 / 6, Self.queryRadius(
      orbital: flerovium.spinNeutralOrbitals[2]), accuracy: 3e-4)
    XCTAssertEqual(2 / 106.828 * 5 / 6, Self.queryRadius(
      orbital: flerovium.spinNeutralOrbitals[3]), accuracy: 3e-4)
    XCTAssertEqual(2 / 106.828 * 5 / 6, Self.queryRadius(
      orbital: flerovium.spinNeutralOrbitals[4]), accuracy: 3e-4)
    XCTAssertEqual(0.033228643, Self.queryRadius(
      orbital: flerovium.spinNeutralOrbitals[5]), accuracy: 1e-3)
    XCTAssertEqual(0.030767119, Self.queryRadius(
      orbital: flerovium.spinNeutralOrbitals[6]), accuracy: 1e-3)
    XCTAssertEqual(0.030767119, Self.queryRadius(
      orbital: flerovium.spinNeutralOrbitals[7]), accuracy: 1e-3)
    XCTAssertEqual(0.030767119, Self.queryRadius(
      orbital: flerovium.spinNeutralOrbitals[8]), accuracy: 1e-3)
    XCTAssertEqual(0.02584201, Self.queryRadius(
      orbital: flerovium.spinNeutralOrbitals[9]), accuracy: 1e-3)
    XCTAssertEqual(0.02584201, Self.queryRadius(
      orbital: flerovium.spinNeutralOrbitals[10]), accuracy: 1e-3)
    XCTAssertEqual(0.025841631, Self.queryRadius(
      orbital: flerovium.spinNeutralOrbitals[11]), accuracy: 1e-3)
    XCTAssertEqual(0.02584201, Self.queryRadius(
      orbital: flerovium.spinNeutralOrbitals[12]), accuracy: 1e-3)
    XCTAssertEqual(0.02584289, Self.queryRadius(
      orbital: flerovium.spinNeutralOrbitals[13]), accuracy: 1e-3)
    XCTAssertEqual(0.049465608, Self.queryRadius(
      orbital: flerovium.spinNeutralOrbitals[14]), accuracy: 1e-3)
    XCTAssertEqual(0.047402453, Self.queryRadius(
      orbital: flerovium.spinNeutralOrbitals[15]), accuracy: 1e-3)
    XCTAssertEqual(0.047402453, Self.queryRadius(
      orbital: flerovium.spinNeutralOrbitals[16]), accuracy: 1e-3)
    XCTAssertEqual(0.047402453, Self.queryRadius(
      orbital: flerovium.spinNeutralOrbitals[17]), accuracy: 1e-3)
  }
}
