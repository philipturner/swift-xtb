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
      octreeFragmentCount += mask64.nonzeroBitCount
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
    descriptor.fragmentCount = 625
    descriptor.position = .zero
    descriptor.sizeExponent = 4
    
    // Proton.
    descriptor.netCharge = +1
    descriptor.netSpin = 0
    let proton = Ansatz(descriptor: descriptor)
    XCTAssertEqual(proton.orbitals.count, 0)
    XCTAssertEqual(proton.occupations, [])
    
    // Neutral hydrogen atom.
    descriptor.netCharge = 0
    descriptor.netSpin = 0.5
    let hydrogenAtom = Ansatz(descriptor: descriptor)
    XCTAssertEqual(hydrogenAtom.orbitals.count, 1)
    XCTAssertEqual(hydrogenAtom.occupations, [1])
    Self.checkFragments(hydrogenAtom.orbitals[0], 625)
    XCTAssertEqual(1, Self.queryRadius(
      orbital: hydrogenAtom.orbitals[0]), accuracy: 0.02)
    
    // Hydride anion (singlet state).
    descriptor.netCharge = -1
    descriptor.netSpin = 0
    let hydrideSinglet = Ansatz(descriptor: descriptor)
    XCTAssertEqual(hydrideSinglet.orbitals.count, 1)
    XCTAssertEqual(hydrideSinglet.occupations, [2])
    Self.checkFragments(hydrideSinglet.orbitals[0], 625)
    XCTAssertEqual(1, Self.queryRadius(
      orbital: hydrideSinglet.orbitals[0]), accuracy: 0.02)
    
    // Hydride anion (triplet state).
    descriptor.atomicNumber = 1
    descriptor.netCharge = -1
    descriptor.netSpin = 1
    let hydrideTriplet = Ansatz(descriptor: descriptor)
    XCTAssertEqual(hydrideTriplet.orbitals.count, 2)
    XCTAssertEqual(hydrideTriplet.occupations, [1, 1])
    Self.checkFragments(hydrideTriplet.orbitals[0], 625)
    Self.checkFragments(hydrideTriplet.orbitals[1], 625)
    XCTAssertEqual(1, Self.queryRadius(
      orbital: hydrideTriplet.orbitals[0]), accuracy: 0.02)
    XCTAssertEqual(2, Self.queryRadius(
      orbital: hydrideTriplet.orbitals[1]), accuracy: 0.08)
  }
  
  func testLithium() throws {
    var descriptor = AnsatzDescriptor()
    descriptor.atomicNumber = 3
    descriptor.fragmentCount = 625
    descriptor.position = .zero
    descriptor.sizeExponent = 4
    
    // Metallic ion.
    descriptor.netCharge = +1
    descriptor.netSpin = 0
    let lithiumIon = Ansatz(descriptor: descriptor)
    XCTAssertEqual(lithiumIon.orbitals.count, 1)
    XCTAssertEqual(lithiumIon.occupations, [2])
    Self.checkFragments(lithiumIon.orbitals[0], 625)
    XCTAssertEqual(0.414, Self.queryRadius(
      orbital: lithiumIon.orbitals[0]), accuracy: 0.02)
    
    // Neutral atom.
    descriptor.netCharge = 0
    descriptor.netSpin = -0.5
    let lithiumNeutral = Ansatz(descriptor: descriptor)
    XCTAssertEqual(lithiumNeutral.orbitals.count, 2)
    XCTAssertEqual(lithiumNeutral.occupations, [2, 1])
    Self.checkFragments(lithiumNeutral.orbitals[0], 625)
    Self.checkFragments(lithiumNeutral.orbitals[1], 625)
    XCTAssertEqual(0.414, Self.queryRadius(
      orbital: lithiumNeutral.orbitals[0]), accuracy: 0.02)
    XCTAssertEqual(2, Self.queryRadius(
      orbital: lithiumNeutral.orbitals[1]), accuracy: 0.08)
    
    // Spin-3/2 atom.
    descriptor.netCharge = 0
    descriptor.netSpin = 1.5
    let lithiumPolarized = Ansatz(descriptor: descriptor)
    XCTAssertEqual(lithiumPolarized.orbitals.count, 3)
    XCTAssertEqual(lithiumPolarized.occupations, [1, 1, 1])
    Self.checkFragments(lithiumPolarized.orbitals[0], 625)
    Self.checkFragments(lithiumPolarized.orbitals[1], 625)
    Self.checkFragments(lithiumPolarized.orbitals[2], 625)
    XCTAssertEqual(0.333, Self.queryRadius(
      orbital: lithiumPolarized.orbitals[0]), accuracy: 0.02)
    XCTAssertEqual(1.414, Self.queryRadius(
      orbital: lithiumPolarized.orbitals[1]), accuracy: 0.03)
    XCTAssertEqual(1.414 * 5 / 6, Self.queryRadius(
      orbital: lithiumPolarized.orbitals[2]), accuracy: 0.03)
  }
  
  // Test chromium.
  // - Test the limitation of GFN-xTB, which doesn't support the true spin
  //   polarization of chromium d-orbitals. Instead, it does 3d^4, 4s^2.
  func testChromium() throws {
    var descriptor = AnsatzDescriptor()
    descriptor.atomicNumber = 24
    descriptor.fragmentCount = 625
    descriptor.position = .zero
    descriptor.sizeExponent = 4
    
    // 3+ oxidation state, with minimal spin polarization.
    descriptor.netCharge = +3
    descriptor.netSpin = 0.5
    let chromiumIon = Ansatz(descriptor: descriptor)
    XCTAssertEqual(chromiumIon.orbitals.count, 11)
    XCTAssertEqual(chromiumIon.occupations, [
      2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2
    ])
    for orbital in chromiumIon.orbitals {
      Self.checkFragments(orbital, 625)
    }
    XCTAssertEqual(1 / 23.414, Self.queryRadius(
      orbital: chromiumIon.orbitals[0]), accuracy: 2e-3)
    XCTAssertEqual(2 / 16.828, Self.queryRadius(
      orbital: chromiumIon.orbitals[1]), accuracy: 4e-3)
    XCTAssertEqual(2 / 16.828 * 5 / 6, Self.queryRadius(
      orbital: chromiumIon.orbitals[2]), accuracy: 4e-3)
    XCTAssertEqual(2 / 16.828 * 5 / 6, Self.queryRadius(
      orbital: chromiumIon.orbitals[3]), accuracy: 4e-3)
    XCTAssertEqual(2 / 16.828 * 5 / 6, Self.queryRadius(
      orbital: chromiumIon.orbitals[4]), accuracy: 4e-3)
    XCTAssertEqual(0.375, Self.queryRadius(
      orbital: chromiumIon.orbitals[5]), accuracy: 0.01)
    XCTAssertEqual(0.375 * 12 / 13, Self.queryRadius(
      orbital: chromiumIon.orbitals[6]), accuracy: 0.01)
    XCTAssertEqual(0.375 * 12 / 13, Self.queryRadius(
      orbital: chromiumIon.orbitals[7]), accuracy: 0.01)
    XCTAssertEqual(0.375 * 12 / 13, Self.queryRadius(
      orbital: chromiumIon.orbitals[8]), accuracy: 0.01)
    XCTAssertEqual(0.2906231, Self.queryRadius(
      orbital: chromiumIon.orbitals[9]), accuracy: 0.01)
    XCTAssertEqual(4 / 4.414, Self.queryRadius(
      orbital: chromiumIon.orbitals[10]), accuracy: 0.04)
    
    
    // Aufbau polarized atom.
    descriptor.netCharge = 0
    descriptor.netSpin = 2
    let chromiumAufbau = Ansatz(descriptor: descriptor)
    XCTAssertEqual(chromiumAufbau.orbitals.count, 14)
    XCTAssertEqual(chromiumAufbau.occupations, [
      2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 2
    ])
    for orbital in chromiumAufbau.orbitals {
      Self.checkFragments(orbital, 625)
    }
    
    // Actual polarized atom.
    descriptor.netCharge = 0
    descriptor.netSpin = 3
    let chromiumActual = Ansatz(descriptor: descriptor)
    XCTAssertEqual(chromiumActual.orbitals.count, 15)
    XCTAssertEqual(chromiumActual.occupations, [
      2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1
    ])
    for orbital in chromiumActual.orbitals {
      Self.checkFragments(orbital, 625)
    }
  }
  
  // Test group (IV) atoms.
  func testGroupIV() throws {
    var descriptor = AnsatzDescriptor()
    descriptor.netCharge = 0
    descriptor.netSpin = 1
    descriptor.position = .zero
    descriptor.sizeExponent = 4
    descriptor.fragmentCount = 625
    
    // Carbon atom (diradical).
    /*
     Converged (1,000,000 fragments)
     expectation radius: 0.18462892
     expectation radius: 0.99920744
     */
    descriptor.atomicNumber = 6
    let carbon = Ansatz(descriptor: descriptor)
    XCTAssertEqual(carbon.orbitals.count, 4)
    XCTAssertEqual(carbon.occupations, [2, 2, 1, 1])
    for orbital in carbon.orbitals {
      Self.checkFragments(orbital, 625)
    }
    XCTAssertEqual(1 / 5.414, Self.queryRadius(
      orbital: carbon.orbitals[0]), accuracy: 0.01)
    XCTAssertEqual(1, Self.queryRadius(
      orbital: carbon.orbitals[1]), accuracy: 0.04)
    
    // Silicon atom (diradical).
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
    XCTAssertEqual(silicon.orbitals.count, 8)
    XCTAssertEqual(silicon.occupations, [2, 2, 2, 2, 2, 2, 1, 1])
    for orbital in silicon.orbitals {
      Self.checkFragments(orbital, 625)
    }
    XCTAssertEqual(1 / 13.414, Self.queryRadius(
      orbital: silicon.orbitals[0]), accuracy: 3e-3)
    XCTAssertEqual(2 / 6.828, Self.queryRadius(
      orbital: silicon.orbitals[1]), accuracy: 6e-3)
    XCTAssertEqual(2 / 6.828 * 5 / 6, Self.queryRadius(
      orbital: silicon.orbitals[2]), accuracy: 6e-3)
    XCTAssertEqual(2 / 6.828 * 5 / 6, Self.queryRadius(
      orbital: silicon.orbitals[3]), accuracy: 6e-3)
    XCTAssertEqual(2 / 6.828 * 5 / 6, Self.queryRadius(
      orbital: silicon.orbitals[4]), accuracy: 6e-3)
    XCTAssertEqual(1.5, Self.queryRadius(
      orbital: silicon.orbitals[5]), accuracy: 0.04)
    
    // Germanium atom (diradical).
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
    XCTAssertEqual(germanium.orbitals.count, 17)
    XCTAssertEqual(germanium.occupations[0], 2)
    XCTAssertEqual(germanium.occupations[14], 2)
    XCTAssertEqual(germanium.occupations[15], 1)
    XCTAssertEqual(germanium.occupations[16], 1)
    for orbital in germanium.orbitals {
      Self.checkFragments(orbital, 625)
    }
    XCTAssertEqual(1 / 31.414, Self.queryRadius(
      orbital: germanium.orbitals[0]), accuracy: 1e-3)
    XCTAssertEqual(2 / 24.828, Self.queryRadius(
      orbital: germanium.orbitals[1]), accuracy: 4e-3)
    XCTAssertEqual(2 / 24.828 * 5 / 6, Self.queryRadius(
      orbital: germanium.orbitals[2]), accuracy: 4e-3)
    XCTAssertEqual(2 / 24.828 * 5 / 6, Self.queryRadius(
      orbital: germanium.orbitals[3]), accuracy: 4e-3)
    XCTAssertEqual(2 / 24.828 * 5 / 6, Self.queryRadius(
      orbital: germanium.orbitals[4]), accuracy: 4e-3)
    XCTAssertEqual(0.36382768, Self.queryRadius(
      orbital: germanium.orbitals[5]), accuracy: 0.02)
    XCTAssertEqual(0.3367374, Self.queryRadius(
      orbital: germanium.orbitals[6]), accuracy: 0.02)
    XCTAssertEqual(0.3367374, Self.queryRadius(
      orbital: germanium.orbitals[7]), accuracy: 0.02)
    XCTAssertEqual(0.3367374, Self.queryRadius(
      orbital: germanium.orbitals[8]), accuracy: 0.02)
    XCTAssertEqual(0.28296816, Self.queryRadius(
      orbital: germanium.orbitals[9]), accuracy: 0.02)
    XCTAssertEqual(0.28296816, Self.queryRadius(
      orbital: germanium.orbitals[10]), accuracy: 0.02)
    XCTAssertEqual(0.2829506, Self.queryRadius(
      orbital: germanium.orbitals[11]), accuracy: 0.02)
    XCTAssertEqual(0.28296816, Self.queryRadius(
      orbital: germanium.orbitals[12]), accuracy: 0.02)
    XCTAssertEqual(0.28294447, Self.queryRadius(
      orbital: germanium.orbitals[13]), accuracy: 0.02)
    XCTAssertEqual(2, Self.queryRadius(
      orbital: germanium.orbitals[14]), accuracy: 0.08)
    
    // Tin atom (diradical).
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
    descriptor.atomicNumber = 50
    descriptor.fragmentCount = 750
    let tin = Ansatz(descriptor: descriptor)
    XCTAssertEqual(tin.orbitals.count, 26)
    XCTAssertEqual(tin.occupations[0], 2)
    XCTAssertEqual(tin.occupations[23], 2)
    XCTAssertEqual(tin.occupations[24], 1)
    XCTAssertEqual(tin.occupations[25], 1)
    for orbital in tin.orbitals {
      Self.checkFragments(orbital, 750)
    }
    XCTAssertEqual(1 / 49.414, Self.queryRadius(
      orbital: tin.orbitals[0]), accuracy: 5e-4)
    XCTAssertEqual(2 / 42.828, Self.queryRadius(
      orbital: tin.orbitals[1]), accuracy: 1e-3)
    XCTAssertEqual(2 / 42.828 * 5 / 6, Self.queryRadius(
      orbital: tin.orbitals[2]), accuracy: 1e-3)
    XCTAssertEqual(2 / 42.828 * 5 / 6, Self.queryRadius(
      orbital: tin.orbitals[3]), accuracy: 1e-3)
    XCTAssertEqual(2 / 42.828 * 5 / 6, Self.queryRadius(
      orbital: tin.orbitals[4]), accuracy: 1e-3)
    XCTAssertEqual(0.114265166, Self.queryRadius(
      orbital: tin.orbitals[5]), accuracy: 6e-3)
    XCTAssertEqual(0.105806686, Self.queryRadius(
      orbital: tin.orbitals[6]), accuracy: 6e-3)
    XCTAssertEqual(0.105806686, Self.queryRadius(
      orbital: tin.orbitals[7]), accuracy: 6e-3)
    XCTAssertEqual(0.105806686, Self.queryRadius(
      orbital: tin.orbitals[8]), accuracy: 6e-3)
    XCTAssertEqual(0.08886481, Self.queryRadius(
      orbital: tin.orbitals[9]), accuracy: 6e-3)
    XCTAssertEqual(0.08886481, Self.queryRadius(
      orbital: tin.orbitals[10]), accuracy: 6e-3)
    XCTAssertEqual(0.08887705, Self.queryRadius(
      orbital: tin.orbitals[11]), accuracy: 6e-3)
    XCTAssertEqual(0.08886481, Self.queryRadius(
      orbital: tin.orbitals[12]), accuracy: 6e-3)
    XCTAssertEqual(0.088864505, Self.queryRadius(
      orbital: tin.orbitals[13]), accuracy: 6e-3)
    XCTAssertEqual(2.5, Self.queryRadius(
      orbital: tin.orbitals[23]), accuracy: 0.06)
    
    // Lead atom (diradical).
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
    descriptor.atomicNumber = 82
    let lead = Ansatz(descriptor: descriptor)
    XCTAssertEqual(lead.orbitals.count, 42)
    XCTAssertEqual(lead.occupations[0], 2)
    XCTAssertEqual(lead.occupations[39], 2)
    XCTAssertEqual(lead.occupations[40], 1)
    XCTAssertEqual(lead.occupations[41], 1)
    for orbital in lead.orbitals {
      Self.checkFragments(orbital, 750)
    }
    XCTAssertEqual(1 / 81.414, Self.queryRadius(
      orbital: lead.orbitals[0]), accuracy: 3e-4)
    XCTAssertEqual(2 / 74.828, Self.queryRadius(
      orbital: lead.orbitals[1]), accuracy: 1e-3)
    XCTAssertEqual(2 / 74.828 * 5 / 6, Self.queryRadius(
      orbital: lead.orbitals[2]), accuracy: 1e-3)
    XCTAssertEqual(2 / 74.828 * 5 / 6, Self.queryRadius(
      orbital: lead.orbitals[3]), accuracy: 1e-3)
    XCTAssertEqual(2 / 74.828 * 5 / 6, Self.queryRadius(
      orbital: lead.orbitals[4]), accuracy: 1e-3)
    XCTAssertEqual(0.051483016, Self.queryRadius(
      orbital: lead.orbitals[5]), accuracy: 3e-3)
    XCTAssertEqual(0.04767065, Self.queryRadius(
      orbital: lead.orbitals[6]), accuracy: 3e-3)
    XCTAssertEqual(0.04767065, Self.queryRadius(
      orbital: lead.orbitals[7]), accuracy: 3e-3)
    XCTAssertEqual(0.04767065, Self.queryRadius(
      orbital: lead.orbitals[8]), accuracy: 3e-3)
    XCTAssertEqual(0.0400402, Self.queryRadius(
      orbital: lead.orbitals[9]), accuracy: 3e-3)
    XCTAssertEqual(0.0400402, Self.queryRadius(
      orbital: lead.orbitals[10]), accuracy: 3e-3)
    XCTAssertEqual(0.040046137, Self.queryRadius(
      orbital: lead.orbitals[11]), accuracy: 3e-3)
    XCTAssertEqual(0.0400402, Self.queryRadius(
      orbital: lead.orbitals[12]), accuracy: 3e-3)
    XCTAssertEqual(0.040042643, Self.queryRadius(
      orbital: lead.orbitals[13]), accuracy: 3e-3)
    XCTAssertEqual(3, Self.queryRadius(
      orbital: lead.orbitals[39]), accuracy: 0.07)
    
    // Flerovium atom (ionized).
    //
    // Stress test for the ansatz code: process the superheavy element
    // flerovium. Its orbitals span many length scales, triggering severe
    // edge cases.
    //
    // We need 5,000 fragments to have good expectation radii.
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
    descriptor.atomicNumber = 114
    descriptor.netCharge = 78
    descriptor.netSpin = 0
    descriptor.fragmentCount = 1500
    let flerovium = Ansatz(descriptor: descriptor)
    XCTAssertEqual(flerovium.orbitals.count, 18)
    XCTAssertEqual(flerovium.occupations[0], 2)
    XCTAssertEqual(flerovium.occupations[17], 2)
    for orbital in flerovium.orbitals {
      Self.checkFragments(orbital, 1500)
    }
    XCTAssertEqual(1 / 113.414, Self.queryRadius(
      orbital: flerovium.orbitals[0]), accuracy: 2e-4)
    XCTAssertEqual(2 / 106.828, Self.queryRadius(
      orbital: flerovium.orbitals[1]), accuracy: 3e-4)
    XCTAssertEqual(2 / 106.828 * 5 / 6, Self.queryRadius(
      orbital: flerovium.orbitals[2]), accuracy: 3e-4)
    XCTAssertEqual(2 / 106.828 * 5 / 6, Self.queryRadius(
      orbital: flerovium.orbitals[3]), accuracy: 3e-4)
    XCTAssertEqual(2 / 106.828 * 5 / 6, Self.queryRadius(
      orbital: flerovium.orbitals[4]), accuracy: 3e-4)
    XCTAssertEqual(0.033228643, Self.queryRadius(
      orbital: flerovium.orbitals[5]), accuracy: 1e-3)
    XCTAssertEqual(0.030767119, Self.queryRadius(
      orbital: flerovium.orbitals[6]), accuracy: 1e-3)
    XCTAssertEqual(0.030767119, Self.queryRadius(
      orbital: flerovium.orbitals[7]), accuracy: 1e-3)
    XCTAssertEqual(0.030767119, Self.queryRadius(
      orbital: flerovium.orbitals[8]), accuracy: 1e-3)
    XCTAssertEqual(0.02584201, Self.queryRadius(
      orbital: flerovium.orbitals[9]), accuracy: 1e-3)
    XCTAssertEqual(0.02584201, Self.queryRadius(
      orbital: flerovium.orbitals[10]), accuracy: 1e-3)
    XCTAssertEqual(0.025841631, Self.queryRadius(
      orbital: flerovium.orbitals[11]), accuracy: 1e-3)
    XCTAssertEqual(0.02584201, Self.queryRadius(
      orbital: flerovium.orbitals[12]), accuracy: 1e-3)
    XCTAssertEqual(0.02584289, Self.queryRadius(
      orbital: flerovium.orbitals[13]), accuracy: 1e-3)
    XCTAssertEqual(0.049465608, Self.queryRadius(
      orbital: flerovium.orbitals[14]), accuracy: 1e-3)
    XCTAssertEqual(0.047402453, Self.queryRadius(
      orbital: flerovium.orbitals[15]), accuracy: 1e-3)
    XCTAssertEqual(0.047402453, Self.queryRadius(
      orbital: flerovium.orbitals[16]), accuracy: 1e-3)
    XCTAssertEqual(0.047402453, Self.queryRadius(
      orbital: flerovium.orbitals[17]), accuracy: 1e-3)
  }
}
