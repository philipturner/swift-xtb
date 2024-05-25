// Test the nitrogen dimer.
func testNitrogen() throws {
  var descriptor = AnsatzDescriptor()
  descriptor.atomicNumbers = [7, 7]
  descriptor.fragmentCount = 5000
  descriptor.positions = [SIMD3(-1, 0, 0), SIMD3(1, 0, 0)]
  descriptor.sizeExponent = 4
  
  // N+ and N-.
  descriptor.netCharges = [+1, -1]
  descriptor.netSpinPolarizations = [0, 0]
  let nitrogenCharged = Ansatz(descriptor: descriptor)
  XCTAssertEqual(nitrogenCharged.spinDownOrbitals.count, 0)
  XCTAssertEqual(nitrogenCharged.spinNeutralOrbitals.count, 7)
  XCTAssertEqual(nitrogenCharged.spinUpOrbitals.count, 0)
  for i in 0..<7 {
    Self.checkFragments(nitrogenCharged.spinNeutralOrbitals[i], 5000)
  }
  XCTAssertEqual(1 / 6.414, Self.queryRadius(
    waveFunction: nitrogenCharged.spinNeutralOrbitals[0],
    nucleusPosition: SIMD3(-1, 0, 0)), accuracy: 5e-3)
  XCTAssertEqual(0.667, Self.queryRadius(
    waveFunction: nitrogenCharged.spinNeutralOrbitals[1],
    nucleusPosition: SIMD3(-1, 0, 0)), accuracy: 0.08)
  XCTAssertEqual(0.667 * 5 / 6, Self.queryRadius(
    waveFunction: nitrogenCharged.spinNeutralOrbitals[2],
    nucleusPosition: SIMD3(-1, 0, 0)), accuracy: 0.08)
  XCTAssertEqual(1 / 6.414, Self.queryRadius(
    waveFunction: nitrogenCharged.spinNeutralOrbitals[3],
    nucleusPosition: SIMD3(1, 0, 0)), accuracy: 5e-3)
  XCTAssertEqual(2 / 1.449, Self.queryRadius(
    waveFunction: nitrogenCharged.spinNeutralOrbitals[4],
    nucleusPosition: SIMD3(1, 0, 0)), accuracy: 0.04)
  XCTAssertEqual(2 / 1.449 * 5 / 6, Self.queryRadius(
    waveFunction: nitrogenCharged.spinNeutralOrbitals[5],
    nucleusPosition: SIMD3(1, 0, 0)), accuracy: 0.04)
  XCTAssertEqual(2 / 1.449 * 5 / 6, Self.queryRadius(
    waveFunction: nitrogenCharged.spinNeutralOrbitals[6],
    nucleusPosition: SIMD3(1, 0, 0)), accuracy: 0.04)
  
  // N(3+) anions.
  // - Literature value for ion atomic radius is 140 pm, or 2.65 Bohr.
  // - The ansatz produces 2.00 Bohr - not bad.
  descriptor.netCharges = [-3, -3]
  descriptor.netSpinPolarizations = [0, 0]
  let nitride = Ansatz(descriptor: descriptor)
  XCTAssertEqual(nitride.spinDownOrbitals.count, 0)
  XCTAssertEqual(nitride.spinNeutralOrbitals.count, 10)
  XCTAssertEqual(nitride.spinUpOrbitals.count, 0)
  for i in 0..<10 {
    Self.checkFragments(nitride.spinNeutralOrbitals[i], 5000)
  }
  XCTAssertEqual(1 / 6.414, Self.queryRadius(
    waveFunction: nitride.spinNeutralOrbitals[0],
    nucleusPosition: SIMD3(-1, 0, 0)), accuracy: 5e-3)
  XCTAssertEqual(2, Self.queryRadius(
    waveFunction: nitride.spinNeutralOrbitals[1],
    nucleusPosition: SIMD3(-1, 0, 0)), accuracy: 0.06)
  XCTAssertEqual(2 * 5 / 6, Self.queryRadius(
    waveFunction: nitride.spinNeutralOrbitals[2],
    nucleusPosition: SIMD3(-1, 0, 0)), accuracy: 0.06)
  XCTAssertEqual(2 * 5 / 6, Self.queryRadius(
    waveFunction: nitride.spinNeutralOrbitals[3],
    nucleusPosition: SIMD3(-1, 0, 0)), accuracy: 0.06)
  XCTAssertEqual(2 * 5 / 6, Self.queryRadius(
    waveFunction: nitride.spinNeutralOrbitals[4],
    nucleusPosition: SIMD3(-1, 0, 0)), accuracy: 0.06)
  XCTAssertEqual(1 / 6.414, Self.queryRadius(
    waveFunction: nitride.spinNeutralOrbitals[5],
    nucleusPosition: SIMD3(1, 0, 0)), accuracy: 5e-3)
  XCTAssertEqual(2, Self.queryRadius(
    waveFunction: nitride.spinNeutralOrbitals[6],
    nucleusPosition: SIMD3(1, 0, 0)), accuracy: 0.06)
  XCTAssertEqual(2 * 5 / 6, Self.queryRadius(
    waveFunction: nitride.spinNeutralOrbitals[7],
    nucleusPosition: SIMD3(1, 0, 0)), accuracy: 0.06)
  XCTAssertEqual(2 * 5 / 6, Self.queryRadius(
    waveFunction: nitride.spinNeutralOrbitals[8],
    nucleusPosition: SIMD3(1, 0, 0)), accuracy: 0.06)
  XCTAssertEqual(2 * 5 / 6, Self.queryRadius(
    waveFunction: nitride.spinNeutralOrbitals[9],
    nucleusPosition: SIMD3(1, 0, 0)), accuracy: 0.06)
  
  // Triplet state.
  descriptor.netCharges = [0, 0]
  descriptor.netSpinPolarizations = [1, 1]
  let nitrogenTriplet = Ansatz(descriptor: descriptor)
  XCTAssertEqual(nitrogenTriplet.spinDownOrbitals.count, 0)
  XCTAssertEqual(nitrogenTriplet.spinNeutralOrbitals.count, 6)
  XCTAssertEqual(nitrogenTriplet.spinUpOrbitals.count, 2)
  for i in 0..<6 {
    Self.checkFragments(nitrogenTriplet.spinNeutralOrbitals[i], 5000)
  }
  Self.checkFragments(nitrogenTriplet.spinUpOrbitals[0], 5000)
  Self.checkFragments(nitrogenTriplet.spinUpOrbitals[1], 5000)
  
  // ROHF simulation that forms 3 covalent bonds.
  descriptor.netCharges = [0, 0]
  descriptor.netSpinPolarizations = [3, -3]
  let nitrogenROHF = Ansatz(descriptor: descriptor)
  XCTAssertEqual(nitrogenROHF.spinDownOrbitals.count, 3)
  XCTAssertEqual(nitrogenROHF.spinNeutralOrbitals.count, 4)
  XCTAssertEqual(nitrogenROHF.spinUpOrbitals.count, 3)
  for i in 0..<4 {
    Self.checkFragments(nitrogenROHF.spinNeutralOrbitals[i], 5000)
  }
  for i in 0..<3 {
    Self.checkFragments(nitrogenROHF.spinDownOrbitals[i], 5000)
    Self.checkFragments(nitrogenROHF.spinUpOrbitals[i], 5000)
  }
}
