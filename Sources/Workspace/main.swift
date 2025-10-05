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
let system: [SIMD4<Float>] = diamondSystem222

// Create the calculator.
var calculatorDesc = xTB_CalculatorDescriptor()
calculatorDesc.atomicNumbers = system.map { UInt8($0.w) }
calculatorDesc.positions = system.map {
  SIMD3($0.x, $0.y, $0.z)
}
calculatorDesc.hamiltonian = .tightBinding
let calculator = xTB_Calculator(descriptor: calculatorDesc)

let energy = calculator.energy
let formattedEnergy = String(format: "%.3f", energy)
print("energy:", formattedEnergy, "zJ")
