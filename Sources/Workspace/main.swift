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

// Load the xTB library with a configuration for maximum performance.
setenv("OMP_STACKSIZE", "2G", 1)
setenv("OMP_NUM_THREADS", "8", 1)
xTB_Library.useLibrary(
  at: "/Users/philipturner/Documents/OpenMM/bypass_dependencies/libxtb.6.dylib")
try! xTB_Library.loadLibrary()
