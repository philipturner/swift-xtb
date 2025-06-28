# Swift Bindings for xTB

## TODO List

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
- Benchmark performance of CPU code in MM4 repo on Windows, ensure no regressions from macOS
  - Use performance test from the test suite, which will be migrated to the new Swift Testing framework
  - Record performance on macOS before and after the migration, to ensure no regressions

Cleanups to all code bases:
- Migrate to Swift 6 (MM4, xTB)
- Migrate tests to the official Swift Testing repo for Swift 6 (MM4)
- Migrate archived code (MM4, xTB) and HardwareCatalog (molecular-renderer) to a dedicated repo
