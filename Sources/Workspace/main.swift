//
//  main.swift
//  xTB
//
//  Created by Philip Turner on 6/19/25.
//

import Foundation
import xTB

// TODO: Make a new script, 'run.sh', that queries the core count and sets
// these environment variables. Then, it calls 'swift run -Xswiftc -Ounchecked'.

let cString1 = getenv("OMP_STACKSIZE")
let cString2 = getenv("OMP_NUM_THREADS")
if let cString1 {
  print("OMP_STACKSIZE:", String(cString: cString1))
}
if let cString2 {
  print("OMP_NUM_THREADS:", String(cString: cString2))
}

// Prepare the environment for maximum performance with xTB.
setenv("OMP_STACKSIZE", "2G", 1)
setenv("OMP_NUM_THREADS", "8", 1) // replace '8' with number of P-cores

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
let system: [SIMD4<Float>] = diamondSystem222

// Create the calculator.
var calculatorDesc = xTB_CalculatorDescriptor()
calculatorDesc.atomicNumbers = system.map { UInt8($0.w) }
calculatorDesc.positions = system.map {
  SIMD3($0.x, $0.y, $0.z)
}
let calculator = xTB_Calculator(descriptor: calculatorDesc)

// Run just one loop iteration.
for _ in 0..<1 {
  calculator.molecule.positions = system.map {
    SIMD3($0.x, $0.y, $0.z)
  }
  
  let checkpoint0 = Date()
  let energy = calculator.energy * xTB_HartreePerZJ
  let checkpoint1 = Date()
  let latency = checkpoint1.timeIntervalSince(checkpoint0)
  
  let formattedLatency = String(format: "%.1f", latency * 1000)
  print()
  print("actual latency:", formattedLatency, "ms")
  
  let formattedEnergy = String(format: "%.7f", energy)
  print("energy:", formattedEnergy, "Ha")
  print("orbitals:", calculator.orbitals.count, calculator.orbitals.eigenvalues.count)
}
