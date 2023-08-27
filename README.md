# Nanomechanical Quantum Chemistry

My experiments with performance-portable DFT algorithms

Goal: Combine a few recent breakthroughs in quantum chemistry. Do this with maximum possible CPU utilization and the simplest possible algorithms.
- Effectively universal XC functional (2021)
  - Removes entire libraries dedicated to handwritten XC heuristics
- Dynamic precision for eigensolvers (2023)
  - Unsure which eigensolvers are applicable to finite difference method. If LOBPCG is, use LOBPCG II.
- Real-space method achieving compute cost parity with plane-wave method (2023)
  - Removes orbital basis sets, drastically simplifying the conceptual complexity
  - Removes the need for FFTs, both an additional library dependency and a bottleneck

Dependencies:
- INQ
- OpenMP
- BLAS
- LAPACK
