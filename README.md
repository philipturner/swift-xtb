# Swift bindings for xTB

## TODO List

GFN2-xTB
- Remove the external charges API, as it's unused
- Properly document how to inject the Accelerate symbolic link
  - Validate that the acceleration works with a fresh install
  - Provide a performance test with the three diamond systems
- Properly document how to set up the OpenMP threads and stack size for optimal performance

GFN-FF
- Avoid fetching `xTB_Results` properties that aren't available to GFN-FF
- Prevent the GFN-FF parameters from spilling to the console
  - https://github.com/grimme-lab/xtb/issues/905
  - May require the full C API for xTB environment verbosity
- Properly document how to fix the GFN-FF crash:
  - Clarify why the crash occurs

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
