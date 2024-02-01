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
- No pseudopotentials
  - Core electrons matter to properly calculate relativistic effects.
  - Pseudopotentials have a non-trivial coupling with the XC functional, complicating testing and trustworthiness of results.
  - Restrict usage to Z <= 36. Use a simple, first-order [relativistic correction](https://www.sciencedirect.com/science/article/abs/pii/S016612800000662X) that only holds for low-Z elements.

## Exchange-Correlation

Exchange-correlation functionals and dispersion corrections should be implemented in separate Swift modules. There should be a plugin-like interface for computing the XC and dispersion terms. This architecture allows different potentials to be applied to different areas of the scene. For example, running a compute-intensive potential on a small fraction of the atoms. One could also choose different potentials for metallic and organic matter.

- [DeepMind 2021 XC functional](https://www.science.org/doi/10.1126/science.abj6511) (2021)
  - Module Name: `DM21`
  - More accurate than the B3LYP functional used for mechanosynthesis research.
  - Provide both the DM21 and DM21mu variants, based on independent reviews of DM21.
  - Facilitate computation of matrix multiplications in a user-specified external library (BLAS, MFA, cuDNN, clBLAST).
- Provide D4 dispersion corrections as a standalone Swift library.
  - Module Name: `D4`

## Finite Differencing Formulas

Variable-resolution orbitals require asymmetric finite differencing. Custom finite difference coefficients were derived for the boundaries between multigrid levels.

<div align="center">

### Symmetric Second-Order

$h^2 f_i'' + \frac{h^4}{12} f_i'''' = f_{i-1} - 2f_i + f_{i+1}$

### Symmetric Fourth-Order

$h^2 f_i'' = \frac{4}{3}(f_{i-1} - 2f_i + f_{i+1}) - \frac{1}{12}(f_{i-2} - 2f_i + f_{i+2})$

### Asymmetric First-Order

$c_- = \frac{1}{h_-}$

$c_+ = \frac{1}{h_+}$

$\frac{h_- + h_+}{2}f_i'' + O(h^2)f_i''' + O(h^3)f_i'''' = c_-f_{i-1} + c_+f_{i+1} - (c_- + c_+)f_i$

### Asymmetric Second-Order

$c_- = \frac{1}{h_-}$

$c_1 = \frac{-2h_-^2 + 8h_+^2}{6h_+^3}$


$c_2 = \frac{h_-^2 - h_+^2}{6h_+^3}$

$(\frac{h_-}{2} + \frac{h_-^2}{6h_+} + \frac{h_+}{3})f_i'' + O(h^3)f_i'''' = c_-f_{i-1} + c_1f_{i+1} + c_2f_{i+2} - (c_- + c_1 + c_2)f_i$

### Doubly Asymmetric Second-Order

TODO: Derive the formula for the case where the +1 sample is fine, but the -1 and +2 samples are coarse.

</div>
