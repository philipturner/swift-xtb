# Density Functional Theory

> This repository is a work in progress.

C++ implementation of performance-portable DFT algorithms, eventually to be GPU-accelerated with Metal and OpenCL.

Goal: Combine a few recent breakthroughs in quantum chemistry. Do this with maximum possible CPU utilization and the simplest possible algorithms.
- [Effectively universal XC functional](https://www.science.org/doi/10.1126/science.abj6511) (2021)
  - Removes entire libraries dedicated to handwritten XC heuristics
- [Dynamic precision for eigensolvers](https://pubs.acs.org/doi/10.1021/acs.jctc.2c00983) (2023)
  - Unsure which eigensolvers are applicable to finite difference method. If LOBPCG is, use LOBPCG II.
- [Real-space method achieving compute cost parity with plane-wave method](https://arxiv.org/pdf/2303.01937.pdf) (2023)
  - Removes orbital basis sets, drastically simplifying the conceptual complexity
  - Removes the need for FFTs, both an additional library dependency and a bottleneck

## Requirements

Operating System:
- macOS or Windows (for now)

Dependencies:
- C++, Clang, Bash
- Pseudopod (from INQ)
- OpenMP
- BLAS
- LAPACK
