//
//  main.swift
//  swift-xtb
//
//  Created by Philip Turner on 6/19/25.
//

import Foundation
import xTB



// GFN-FF
//
// System | Correct Energy  |
// ------ | --------------- |
// 122    |   -57126.394 zJ |
// 222    |  -100762.569 zJ |
// 233    |  -208433.751 zJ |
//
// GFN2-xTB
//
// System | Correct Energy  |
// ------ | --------------- |
// 122    |  -452321.592 zJ |
// 222    |  -841344.658 zJ |
// 233    | -1796496.861 zJ |



let cString1 = getenv("OMP_STACKSIZE")
let cString2 = getenv("OMP_NUM_THREADS")
if let cString1 {
  print("OMP_STACKSIZE:", String(cString: cString1))
}
if let cString2 {
  print("OMP_NUM_THREADS:", String(cString: cString2))
}

// Suppress unwanted output from GFN-FF.
let url = FileManager.default.temporaryDirectory
let path = url.relativePath
let worked = FileManager.default.changeCurrentDirectoryPath(path)
guard worked else {
  fatalError("Could not redirect gfnff_topo directory.")
}
xTB_Environment.verbosity = .muted
xTB_Environment.setOutput("/dev/null")

// Select the system.
let system: [SIMD4<Float>] = diamondSystem233

// Create the calculator.
var calculatorDesc = xTB_CalculatorDescriptor()
calculatorDesc.atomicNumbers = system.map { UInt8($0.w) }
calculatorDesc.positions = system.map {
  SIMD3($0.x, $0.y, $0.z)
}
//calculatorDesc.hamiltonian = .forceField
let calculator = xTB_Calculator(descriptor: calculatorDesc)

// Run just one loop iteration.
var minLatency: Double = 1_000_000
for _ in 0..<10 {
  calculator.molecule.positions = system.map {
    SIMD3($0.x, $0.y, $0.z)
  }
  
  let checkpoint0 = Date()
  let energy = calculator.energy
  let checkpoint1 = Date()
  let latency = checkpoint1.timeIntervalSince(checkpoint0)
  minLatency = min(latency, minLatency)
  
  let formattedLatency = String(format: "%.1f", latency * 1000)
  print()
  print("actual latency:", formattedLatency, "ms")
  
  let formattedEnergy = String(format: "%.3f", energy)
  print("energy:", formattedEnergy, "zJ")
  print("orbitals:", calculator.orbitals.count, calculator.orbitals.eigenvalues.count)
  print("bond orders:", calculator.molecule.bondOrders.count)
  print("charges:", calculator.molecule.charges.count)
  print("forces:", calculator.molecule.forces.count)
  
  // If this is commented out, then the energies change.
  print()
  print(xTB_Environment.status)
  print(xTB_Environment.flushErrorStack())
  print(xTB_Environment.status)
  xTB_Environment.show()
  print(xTB_Environment.status)
  print()
}

// Summarize the results
let formattedLatency = String(format: "%.1f", minLatency * 1000)
print()
print("minimum latency:", formattedLatency, "ms")

