import XCTest
import Mechanosynthesis
import Numerics

// Currently a temporary class to facilitate an experiment.
final class GridTests: XCTestCase {
  // Set up the water molecule.
  // - H-O bond length: 95.84 pm
  // - H-O-H angle: 104.45°
  static func createWaterAnsatz() -> Ansatz {
    // Create the oxygen atom.
    var positions: [SIMD3<Float>] = []
    let oxygenPosition: SIMD3<Float> = .zero
    positions.append(oxygenPosition)
    
    // Create the hydrogen atoms.
    for hydrogenID in 0..<2 {
      let angleDegrees = Float(hydrogenID) * 104.45
      let angleRadians = angleDegrees * (Float.pi / 180)
      
      // Create a direction vector from the angle.
      let directionX = Float.cos(angleRadians)
      let directionZ = -Float.sin(angleRadians)
      let direction = SIMD3(directionX, 0, directionZ)
      
      // Determine the hydrogen atom's position.
      let bohrRadiusInPm: Float = 52.9177210903
      let bondLengthInBohr = 95.84 / bohrRadiusInPm
      let hydrogenPosition = oxygenPosition + direction * bondLengthInBohr
      
      // Append the hydrogen to the array.
      positions.append(hydrogenPosition)
    }
    
    // Specify the topology.
    //
    // The nuclear coordinates range from ±1.8 Bohr.
    // The octree spans from ±0.5 * 16.0 Bohr.
    var ansatzDesc = AnsatzDescriptor()
    ansatzDesc.atomicNumbers = [8, 1, 1]
    ansatzDesc.positions = positions
    ansatzDesc.sizeExponent = 4
    
    // Set the quality to ≥1000 fragments/electron.
    ansatzDesc.fragmentCount = 1000
    
    // Specify the net charges and spins.
    ansatzDesc.netCharges = [0, 0, 0]
    ansatzDesc.netSpinPolarizations = [0, 1, -1]
    return Ansatz(descriptor: ansatzDesc)
  }
  
  // Some experimentation with 'Grid', figuring out how to transfer data from
  // the ansatz's 'Octree' to the grid. Will most likely be archived or deleted.
  func testWorkspace() throws {
    let ansatz = Self.createWaterAnsatz()
    
    // Objective: transfer each wavefunction to a dedicated 'Grid'. Then,
    // represent the charge density of the system.
    // - Create a merged grid that is a union of the wavefunctions' grids.
    // - Expose the AtomicOrbital API and query each wavefunction's value on
    //   the higher-resolution merged grid.
    // - This action is "de-localizing" the wavefunctions, forcing them to use
    //   the same dense basis.
    
    // First: find the bounding box of one wavefunction at the 2x2x2 Bohr
    // granularity.
    let waveFunction = ansatz.spinNeutralWaveFunctions[0]
    var boundsMinimum = SIMD3<Float>(repeating: .greatestFiniteMagnitude)
    var boundsMaximum = SIMD3<Float>(repeating: -.greatestFiniteMagnitude)
    for node in waveFunction.octree.nodes {
      guard node.spacing <= 2 else {
        continue
      }
      let nodeMinimum = node.center - node.spacing / 2
      let nodeMaximum = node.center + node.spacing / 2
      boundsMinimum.replace(
        with: nodeMinimum, where: nodeMinimum .< boundsMinimum)
      boundsMaximum.replace(
        with: nodeMaximum, where: nodeMaximum .> boundsMaximum)
    }
    print(boundsMinimum)
    print(boundsMaximum)
    
    var gridDesc = GridDescriptor()
    gridDesc.offset = SIMD3<Int>(boundsMinimum)
    gridDesc.dimensions = SIMD3<Int>(boundsMaximum - boundsMinimum)
    let grid = Grid(descriptor: gridDesc)
    print(grid.offset)
    print(grid.dimensions)
  }
}
