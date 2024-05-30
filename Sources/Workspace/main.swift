//
//  main.swift
//
//
//  Created by Philip Turner on 5/25/24.
//

import Darwin
import xTB

// From the command line:
// swift run -Xswiftc -Ounchecked --debugger Workspace -o run

// Configure xTB for maximum performance.
setenv("OMP_STACKSIZE", "2G", 1)
setenv("OMP_NUM_THREADS", "8", 1)
xTBLibrary.useLibrary(
  at: "/Users/philipturner/Documents/OpenMM/bypass_dependencies/libxtb.6.dylib")
try! xTBLibrary.loadLibrary()

// Check that xTB is loading correctly.
let xtb_getAPIVersion: @convention(c) () -> Int32 =
xTBLibrary.loadSymbol(name: "xtb_getAPIVersion")
print(xtb_getAPIVersion())
