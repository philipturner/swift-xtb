# Density Functional Theory

Goal: Combine a few recent advances in quantum chemistry. Do this with maximum possible CPU utilization and the simplest possible algorithms. Due to its simplicity, all of the code can be ported to OpenCL and Metal.

- Real-space formalism
  - Removes orbital basis sets, drastically simplifying the functional form.
  - Removes FFTs, a bottleneck and library dependency.
  - Most DFT libraries (Gaussian, GAMESS, TeraChem) use the plane-wave formalism. This formalism is well-suited to CPUs, but not GPUs.
- Variable-resolution orbitals to accelerate the onset of $O(n)$ scaling.
  - Loosely constrain each orbital fragment to have the same probability.
- [Dynamic precision for eigensolvers](https://pubs.acs.org/doi/10.1021/acs.jctc.2c00983) (2023)
  - Allows DFT to run on consumer hardware with few FP64 units.
  - Remove LOBPCG and all linear algebra, such as `dsyevd`. Solve the eigenproblem with a linear-scaling algorithm.
- [DeepMind 2021 XC functional](https://www.science.org/doi/10.1126/science.abj6511) (2021)
  - More accurate than the B3LYP functional used for mechanosynthesis research.
  - Provide both the DM21 and DM21mu variants, based on independent reviews of DM21.
  - Provide both D3 and D4 dispersion corrections as the only external dependencies.
- No pseudopotentials
  - Core electrons matter to properly calculate relativistic effects.
  - Pseudopotentials have a non-trivial coupling with the XC functional, complicating testing and trustworthiness of results.
  - Restrict usage to Z <= 36. Use a simple, first-order [relativistic correction](https://www.sciencedirect.com/science/article/abs/pii/S016612800000662X) that only holds for low-Z elements.
