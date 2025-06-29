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
}

// Summarize the results
let formattedLatency = String(format: "%.2f", minLatency * 1000)
print()
print("minimum latency:", formattedLatency, "ms")

