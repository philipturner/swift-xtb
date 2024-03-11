# Dependencies

This page lists current and planned dependencies for the library. It is not exhaustive at the moment.

Dependencies:
- DFT-D4
  - PythonKit-style linking to avoid compiler issues.
- DM21
  - Weights embedded into source tree.
- libxc
  - Linked at compile-time for now.
- OpenCL
  - SIMD-scoped shuffles through either Khronos extensions or assembly injection.
