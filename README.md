# Density Functional Theory

Overview:
- Real-space formalism
  - Removes orbital basis sets, drastically simplifying the functional form.
  - Removes FFTs, a bottleneck and library dependency.
  - Most DFT libraries (Gaussian, GAMESS, TeraChem) use the plane-wave formalism. This formalism is well-suited to CPUs, but not GPUs.
- Variable-resolution orbitals to accelerate the onset of $O(n)$ scaling
  - Loosely constrain each orbital fragment to have the same probability.
  - User specifies the probability density for orbital fragments at the 2x2x2 granularity, where the most expensive operations are performed.
- [Dynamic precision for eigensolvers](https://pubs.acs.org/doi/10.1021/acs.jctc.2c00983) (2023)
  - Allows DFT to run on consumer hardware with few FP64 units.
  - Remove LOBPCG and all linear algebra, such as `dsyevd`. Solve the eigenproblem with a linear-scaling algorithm.
  - Perform all computations in FP32, with compensated summation when necessary.
- No pseudopotentials
  - Core electrons matter to properly calculate relativistic effects.
  - Fixed pseudopotentials have a non-trivial coupling with the XC functional, complicating testing and trustworthiness of results.
  - Generate numerical pseudopotentials at runtime, let the user recycle them for similar computations.
- [DeepMind 2021 XC functional](https://www.science.org/doi/10.1126/science.abj6511) (2021)
  - More accurate than the B3LYP functional used for mechanosynthesis research.
  - Provide the following XC functionals to the user: exchange-only LDA, Hartree-Fock, DM21, DM21mu.
  - Allow the dispersion component of correlation to be disabled, especially if the dylib can't be located.

Dependencies:
- DFT-D4
  - PythonKit-style linking to avoid compiler issues.
- DM21
  - Weights embedded into source tree.
- Metal (only on Apple platforms)
  - MFA binary embedded into source tree, potentially with fused activations.
- OpenCL
  - SIMD-scoped shuffles through either Khronos extensions or assembly injection.

## Finite Differencing

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

$c_1 = \frac{4h_+^2 - h_-^2}{3h_+^3}$

$c_2 = \frac{h_-^2 - h_+^2}{6h_+^3}$

$(\frac{h_-}{2} + \frac{h_-^2 + 2h_+^2}{6h_+})f_i'' + O(h^3)f_i'''' = c_-f_{i-1} + c_1f_{i+1} + c_2f_{i+2} - (c_- + c_1 + c_2)f_i$

$O(h^3) = \frac{1}{12}(\frac{h_-^3}{2} + \frac{ h_+^2 (7h_-^2 - 4h_+^2) }{6h_+})$

### Doubly Asymmetric Second-Order

$c_- = \frac{1}{h_-}$

$c_1 = \frac{(h_1 + h_2)^2 - h_-^2}{h_1 (2h_1h_2 + h_2^2)}$

$c_2 = \frac{h_-^2 - h_1^2}{(h_1 + h_2)(2h_1h_2 + h_2^2)}$

$(\frac{h_-}{2} + \frac{h_-^2 + h_1^2 + h_1h_2}{4h_1 + 2h_2})f_i'' + O(h^3)f_i'''' = c_-f_{i-1} + c_1f_{i+1} + c_2f_{i+2} - (c_- + c_1 + c_2)f_i$

$O(h^3) = \frac{1}{12}(\frac{h_-^3}{2} + \frac{h_-^2 (3h_1^2 + 3h_1h_2 + h_2^2) - h_1^2 (h_1 + h_2)^2 }{4h_1 + 2h_2})$

</div>

## Relativity Corrections

Relativity corrections can only be computed in momentum space (plane-wave formalism). To be more efficient and sparse, this DFT library uses the real-space formalism. The best solution is a [rough approximation](https://doi.org/10.1088/1361-6404/ac0ecc).

<div align="center">

$\left[-\frac{1}{2} \nabla^2 + V\right] \Psi = E \Psi$

$H \Psi = E \Psi $

</div>

The non-relativistic Schrodinger equation (above) transforms into the expression below.

<div align="center">

$\gamma = \sqrt{1 + \langle \Psi | \hat{p}^2 | \Psi \rangle / c^2}$

$\left[-\frac{1}{\gamma + 1} \nabla^2 + V\right] \Psi = E \Psi$

$H \Psi = E \Psi $

</div>

> TODO: Compare the energies and contracted radii to the [exact solution](https://doi.org/10.1038/s41598-020-71505-w).

## Screened Exchange

DM21 is a range-separated exchange-correlation functional. The regular Fock exchange integral can be evaluated with a Poisson solver, but the range-separated term is more challenging. It cannot be directly solved in real-space efficiently.

<div align="center">

$e^{HF}(r) = -\frac{1}{2} \Sigma_i\Sigma_j\Psi_i(r)\Psi_j(r) \int\Psi_i(r')\Psi_j(r') \frac{1}{|r-r'|} d^{3}r'$

$v^{HF}(r) = \int\Psi_i(r')\Psi_j(r') \frac{1}{|r-r'|} d^{3}r'$

$\nabla^2 v^{HF}(r) = -4\pi \Psi_i(r)\Psi_j(r)$

</div>

A simple alternative is replacing the `erf`-type screening with `exp`-type screening. This version can be solved exactly with multigrids. The screening coefficient ($\omega$) is 0.4 Bohr<sup>-1</sup>.

<div align="center">

$v^{\omega HF}(r) = \int\Psi_i(r')\Psi_j(r') \frac{erf(\omega|r-r'|)}{|r-r'|} d^{3}r'$

$v^{\omega HF}(r) \approx \int\Psi_i(r')\Psi_j(r') \frac{\exp(-\omega|r-r'|)}{|r-r'|} d^{3}r'$

$(\nabla^2 - \omega^2) v^{\omega HF}(r) = -4\pi \Psi_i(r)\Psi_j(r)$

</div>

To improve accuracy, the screened function is approximated by a linear combination of three `exp`-type functions. The first is recycled from the regular exchange integral. The third is comparatively cheap. The [overall function](https://www.desmos.com/calculator/oyasyidilk) deviates from `erf` by less than 5%.

<div align="center">

$\frac{erf(\omega x)}{x} \approx \frac{1 - 1.4\exp(-2\omega x) - 0.4\exp(-4.3\omega x)}{x}$

| Component | Poisson Equation | Cutoff |
| :---: | :---: | :---: |
| $\frac{1}{x}$ | $\nabla^2 v(r) = -4\pi \rho(r)$ | none |
| $\frac{\exp(-2\omega x)}{x}$ | $(\nabla^2 - (2\omega)^2) v(r) = -4\pi \rho(r)$ | 4.3 Bohr |
| $\frac{\exp(-4.3\omega x)}{x}$ | $(\nabla^2 - (4.3\omega)^2) v(r) = -4\pi \rho(r)$ | 2.15 Bohr |

</div>

To compute the `erf`-type potential, first solve each component function in an independent Poisson equation. Then, sum the components at each specific grid point. Integral computations beyond a cutoff radius are skipped to improve efficiency. The cutoffs make the function piecewise.
