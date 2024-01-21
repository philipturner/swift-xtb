# Density Functional Theory

Short-Term Goal: Wrap the [GOSPEL](https://gitlab.com/jhwoo15/gospel) Python package in an ergonomic Swift API.

Long-Term Goal: Combine a few recent advances in quantum chemistry. Do this with maximum possible CPU utilization and the simplest possible algorithms. Then, port the most compute-intensive parts to OpenCL.

- Real-space formalism
  - Removes orbital basis sets, drastically simplifying the functional form.
  - Reduces the number of FFTs, a bottleneck that dominates computation time.
  - Most DFT libraries (Gaussian, GAMESS, TeraChem) uses the plane-wave formalism. This formalism is well-suited to CPUs, but not GPUs.
- [Universal XC functional](https://www.science.org/doi/10.1126/science.abj6511) (2021)
  - More accurate than the B3LYP functional used for mechanosynthesis research, or at least not significantly worse.
  - The XC functional is often 90% of the maintenance and complexity of a DFT codebase. DeepMind's neural network makes the XC code ridiculously simple.
- [Dynamic precision for eigensolvers](https://pubs.acs.org/doi/10.1021/acs.jctc.2c00983) (2023)
  - Allows DFT to run on consumer hardware with few FP64 units.
  - Use a similar solver described there, except replacing LOBPCG with LOBPCG II. This reduces bottlenecks from eigendecomposition (`eigh`) by 27x.
