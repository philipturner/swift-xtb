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

TODO: Change the meaning of fragment count.
TODO: Multithread the ansatz generation.
TODO: Change the workspace to a Swift executable module.
