//
//  main.swift
//
//
//  Created by Philip Turner on 5/25/24.
//

import Accelerate
import Foundation
import LinearAlgebra
import Numerics
import xTB

// Copy the dylib from Homebrew Cellar to the folder for hacked dylibs. Use
// 'otool' to replace the OpenBLAS dependency with Accelerate. To do this:
// - copy libxtb.dylib into custom folder
// - otool -L "path to libxtb.dylib"
// - find the address of the OpenBLAS in the output
// - install_name_tool -change "path to libopenblas.dylib" \
//   "/System/Library/Frameworks/Accelerate.framework/Versions/A/Accelerate" \
//   "path to libxtb.dylib"
// Prepare the environment for maximum performance with xTB.
setenv("OMP_STACKSIZE", "2G", 1)
setenv("OMP_NUM_THREADS", "8", 1) // replace '8' with number of P-cores

// Load the 'xtb' dylib.
let pathPart1 = "/Users/philipturner/Documents/OpenMM"
let pathPart2 = "/bypass_dependencies/libxtb.6.dylib"
xTB_Library.useLibrary(at: pathPart1 + pathPart2)
try! xTB_Library.loadLibrary()

xTB_Environment.verbosity = .minimal

var calculatorDesc = xTB_CalculatorDescriptor()
calculatorDesc.atomicNumbers = diamondSystem122.map { UInt8($0.w) }
let calculator = xTB_Calculator(descriptor: calculatorDesc)

for _ in 0..<10 {
  calculator.molecule.positions = diamondSystem122.map {
    SIMD3($0.x, $0.y, $0.z)
  }
  
  let checkpoint0 = Date()
  _ = calculator.energy
  let checkpoint1 = Date()
  let latency = checkpoint1.timeIntervalSince(checkpoint0)
  print()
  print("actual latency:", latency)
  
  guard calculator.orbitals.count == 196 else {
    fatalError("Unexpected problem size: \(calculator.orbitals.count)")
  }
}
