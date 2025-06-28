# Swift Bindings for xTB

## TODO List

GFN2-xTB
- Why are there text files for GFN2-xTB parameters?
  - Try deleting the Homebrew installation and the build artifacts from `compile-xtb.sh`.
  - This may give hints to the GFN-FF crash relating to the `FileManager` directory in Xcode.

GFN-FF
- Properly document how to fix the GFN-FF crash:
  - Clarify why the crash occurs

API improvements:
- Automatically suppress `gfnff_topo` file writing in a robust manner
  - Figure out exactly when it's written, then switch back to the previous directory afterward
  - Purge `gfnff_topo` and `gfnff_charges` from the NSTemporaryDirectory, so that every initialization of `xTB_Calculator` regenerates the GFN-FF parameters from scratch.
  - Make sure to use a cross-platform equivalent of NSTemporaryDirectory
    - Get the xTB bindings working on Windows

New goals:
- Remove the energy minimizer from MM4
- Get the code correctly compiling from source on Windows
  - Benchmark the diamond systems to test for correct optimization flags

Cleanups to all code bases:
- Migrate to Swift 6 (MM4, xTB)
- Migrate tests to the official Swift Testing repo for Swift 6 (MM4)
- Migrate archived code (MM4, xTB) and HardwareCatalog (molecular-renderer) to a dedicated repo

## Current Documentation

```swift

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

// Fix the GFN-FF crash.
FileManager.default.changeCurrentDirectoryPath("/Users/philipturner")

// Load the 'xtb' dylib.
let pathPart1 = "/Users/philipturner/Documents/OpenMM"
let pathPart2 = "/bypass_dependencies/libxtb.6.dylib"
xTB_Library.useLibrary(at: pathPart1 + pathPart2)
try! xTB_Library.loadLibrary()

// Mute the output to the console.
xTB_Environment.verbosity = .muted

```
