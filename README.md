# Swift Bindings for xTB

## TODO List

New goals:
- Remove the energy minimizer from MM4
- Get the code correctly compiling from source on Windows
  - Benchmark the diamond systems to test for correct optimization flags
- Benchmark performance of CPU code in MM4 repo on Windows, ensure no regressions from macOS
  - Use performance test from the test suite, which will be migrated to the new Swift Testing framework
  - Record performance on macOS before and after the migration, to ensure no regressions
- Integrate all of this maintenance into MM4 in a single branch, '2025-cleanups'.

Cleanups to all code bases:
- Migrate to Swift 6 (MM4)
- Migrate tests to the official Swift Testing repo for Swift 6 (MM4)
- Migrate archived code (MM4, xTB) and HardwareCatalog (molecular-renderer) to a dedicated repo
