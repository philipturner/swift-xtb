# Mechanosynthesis

> This library is a work in progress.

Quantum mechanics simulator for molecular nanotechnology.

Goals:
- Designed for GPU acceleration with OpenCL
- Real-space formalism
  - Removes FFTs, an $O(n^2\log{n})$ bottleneck.
  - Most DFT libraries (Gaussian, GAMESS, TeraChem) use the plane-wave/basis-set formalism. This formalism is well-suited to CPUs, but not GPUs.
  - Adaptive mesh refinement at the per-electron granularity, instead of the typical global granularity. This reduces the bottleneck from $O(n^2)$ to $O(n\log{n})$
- [Dynamic precision for eigensolvers](https://pubs.acs.org/doi/10.1021/acs.jctc.2c00983) (2023)
  - Allows DFT to run on consumer hardware with few FP64 units.
  - Perform all computations in FP32, with compensated summation when necessary.
- No pseudopotentials
  - Core electrons matter to properly calculate relativistic effects.
  - Fixed pseudopotentials have a non-trivial coupling with the XC functional, complicating testing and trustworthiness of results.
- [DeepMind 2021 XC functional](https://www.science.org/doi/10.1126/science.abj6511) (2021)
  - More accurate than the B3LYP functional used for mechanosynthesis research.
  - Provide all DM21 variants, in addition to the functionals from libxc.

Documentation:
- [Dependencies](./Documentation/Pages/Dependencies.md)
- [Relativity Corrections](./Documentation/Pages/RelativityCorrections)
- [Screened Exchange](./Documentation/Pages/ScreenedExchange.md)