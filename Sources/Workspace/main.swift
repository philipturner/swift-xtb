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
