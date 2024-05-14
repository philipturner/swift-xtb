# Mechanosynthesis

> This library is a work in progress.

Quantum mechanics simulator for molecular nanotechnology.

Goals:
- Designed for GPU acceleration with OpenCL
  - Perform all computations in FP32, with compensated summation when necessary.
- Real-space formalism
  - Removes FFTs, an $O(n^2\log{n})$ bottleneck.
  - Most DFT libraries (Gaussian, GAMESS, TeraChem) use the plane-wave/basis-set formalism. Basis sets add unnecessary complexity and obfuscate pathways to linear scaling.
- Linear scaling
  - Adaptive mesh refinement at the per-electron granularity, instead of the typical global granularity. This reduces the bottleneck from $O(n^2)$ to $O(n\log{n})$.
  - Linear scaling likely possible for insulators, with localized electrons.
- No pseudopotentials
  - Large-core pseudopotentials interact with the XC functional, polluting simulation results.
  - Similar to the issue with the AO basis: BSSE pollutes simulation results.

Dependencies:
- DFT-D4
  - PythonKit-style linking to avoid compiler issues.
- libxc
  - Linked at compile-time for now.
- OpenCL
  - SIMD-scoped shuffles through either Khronos extensions or assembly injection.

## TODO List for this Branch

libxc
- LDA
  - LDA exchange energy density and potential
  - LDA correlation energy density and potential
  - LSDA
- GGA
  - PBE
  - r<sup>2</sup>SCAN

Electrostatics
- Self-repulsion integral ✅
- Neumann boundaries ✅
- Fuzzy cell decomposition
- Laplacian preconditioner

FiniteDifferencing
- Quadratic interpolation ✅
- Mehrstellen operator
  - Generation of a residual vector
  - Rayleigh quotient estimation
- Coarsening operator

LinearSolver
- Jacobi
- Gauss-Seidel
- Conjugate gradient
- Multigrid

Relativity
- Recursive definition of kinetic energy ✅
- Solving PGPP equation for radial wavefunction of hydrogen-like atoms ✅
