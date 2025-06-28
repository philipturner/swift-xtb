# Swift Bindings for xTB

## TODO List

GFN2-xTB
- Properly document how to inject the Accelerate symbolic link
  - Validate that the acceleration works with a fresh install
  - Provide a performance test with the three diamond systems
  - Unable to quantify performance with `diamondSystem233` due to a crash with the current xTB version
- Properly document how to set up the OpenMP threads and stack size for optimal performance

GFN-FF
- Properly document how to fix the GFN-FF crash:
  - Clarify why the crash occurs
  - Try to reproduce the crash in the old molecular-renderer, with the new xTB bindings
    - Patch up the current main branch of molecular-renderer in a non-main branch. This can serve as a surrogate for testing simulators until the overhaul is complete.
    - Fresh `bypass_dependencies` folder in the fresh MolecularRenderer directory. No need for the extra dylibs (probably due to swift-gif and other community packages). Using the most recent OpenMM and the standard OpenCL backend.

API improvements:
- Update the API for Swift 6
- Automatically suppress `gfnff_topo` file writing in a robust manner
  - Figure out exactly when it's written, then switch back to the previous directory afterward
  - Purge `gfnff_topo` and `gfnff_charges` from the NSTemporaryDirectory, so that every initialization of `xTB_Calculator` regenerates the GFN-FF parameters from scratch.
  - Make sure to use a cross-platform equivalent of NSTemporaryDirectory
    - Get the xTB bindings working on Windows
- Automatically set `OMP_STACKSIZE` and `OMP_NUM_THREADS` prior to invoking either GFN2-xTB or GFN-FF, in a testable manner.
  - Set the number of CPU cores to `perflevel0.physicalcpu` during the `run.sh` script, but only on macOS.
- Ensure all issues currently on the README are addressed. Then, proceed with intercepting linear algebra library calls.

End goals:
- Production-ready API with an opt-in FP32 mode, on both macOS and Windows
- Able to gather data about contributions to latency across a diverse set of environments, for the 3 diamond systems

New goals (just spilling my TODO list):
- Remove the energy minimizer from MM4
- Fix the errors with HDL
- Fix up the xTB bindings and move on
  - Get automatic suppression of `gfnff_topo` working 
  - Get the code correctly compiling from source on Windows, and figure out the need for `param_gfn2.txt`
  - Get benchmarks on Windows
- Create simulators TODO list in molecular-renderer
- Get all simulators integrated into the new molecular-renderer

The latest set of goals combines with all of the concerns preceding it. Resolve all of them during this round of software maintenance.

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
