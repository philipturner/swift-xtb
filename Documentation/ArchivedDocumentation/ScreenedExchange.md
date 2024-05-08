# Screened Exchange

[DeepMind 2021 XC functional](https://www.science.org/doi/10.1126/science.abj6511) (2021)
- More accurate than the B3LYP functional used for mechanosynthesis research.
- Provide all DM21 variants, in addition to the functionals from libxc.

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

$(\nabla^2 - \omega^2) v^{\omega HF}(r) \approx -4\pi \Psi_i(r)\Psi_j(r)$

</div>

To improve accuracy, the screened function is approximated by a linear combination of three `exp`-type functions. The first is recycled from the regular exchange integral. The third is comparatively cheap. The [overall function](https://www.desmos.com/calculator/ukk2k5z816) deviates from `erf` by less than 5%.

<div align="center">

$\frac{erf(\omega x)}{x} \approx \frac{1 - 1.4\exp(-2\omega x) + 0.4\exp(-4.3\omega x)}{x}$

| Component | Poisson Equation | Cutoff |
| :---: | :---: | :---: |
| $\frac{1}{x}$ | $\nabla^2 v(r) = -4\pi \rho(r)$ | none |
| $\frac{\exp(-2\omega x)}{x}$ | $(\nabla^2 - (2\omega)^2) v(r) = -4\pi \rho(r)$ | 4.3 Bohr |
| $\frac{\exp(-4.3\omega x)}{x}$ | $(\nabla^2 - (4.3\omega)^2) v(r) = -4\pi \rho(r)$ | 2.15 Bohr |

</div>

To compute the `erf`-type potential, first solve each component function in an independent Poisson equation. Then, sum the components at each specific grid point. Integral computations beyond a cutoff radius are skipped to improve efficiency. The cutoffs make the function piecewise.

