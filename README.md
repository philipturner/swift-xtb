# Density Functional Theory

Overview:

- Real-space formalism
  - Removes orbital basis sets, drastically simplifying the functional form.
  - Removes FFTs, a bottleneck and library dependency.
  - Most DFT libraries (Gaussian, GAMESS, TeraChem) use the plane-wave formalism. This formalism is well-suited to CPUs, but not GPUs.
- Variable-resolution orbitals to accelerate the onset of $O(n)$ scaling
  - Loosely constrain each orbital fragment to have the same probability.
- [Dynamic precision for eigensolvers](https://pubs.acs.org/doi/10.1021/acs.jctc.2c00983) (2023)
  - Allows DFT to run on consumer hardware with few FP64 units.
  - Remove LOBPCG and all linear algebra, such as `dsyevd`. Solve the eigenproblem with a linear-scaling algorithm.
  - Perform 100% of computations in FP32, with compensated summation when necessary.
- No pseudopotentials
  - Fixed pseudopotentials have a non-trivial coupling with the XC functional, complicating testing and trustworthiness of results. Find a general-purpose alternative that generates pseudopotentials at runtime.
  - Core electrons matter to properly calculate relativistic effects.
- No external dependencies except OpenCL
  - Requires conformance to the OpenCL 2 extension for sub-group shuffles and reductions.
  - Apple silicon conforms through [AIR workaround](https://github.com/philipturner/opencl-metal-stdlib).
  - Nvidia might require injecting PTX assembly.

## Exchange-Correlation

Exchange-correlation functionals and dispersion corrections should be implemented in separate Swift modules. There should be a programmable interface for computing the XC and dispersion terms. The workflow should not require a dependency on LibXC.

- [DeepMind 2021 XC functional](https://www.science.org/doi/10.1126/science.abj6511) (2021)
  - Module Name: `DM21`
  - More accurate than the B3LYP functional used for mechanosynthesis research.
  - Provide all four DM21 variants.
  - Run matrix multiplications in vendor-specific libraries (Accelerate, MFA, cuBLAS, clBLAST, etc.).
- Provide D4 dispersion corrections as a standalone Swift library.
  - Module Name: `D4`
  - Create a Bash script that downloads a Fortran compiler and builds the library from source.
- The `DFT` module encapsulates the primitives for XC functionals.
  - Density in each spin channel.
  - First and second derivatives of density in each spin channel.
  - Kinetic energy density in each spin channel.
  - Exact exchange, with a range separation parameter. This may be requested multiple times with differing parameters.

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

$H \Psi = E \Psi $

$\left[-\frac{1}{2} \nabla^2 + V\right] \Psi = E \Psi$

</div>

The non-relativistic Schrodinger equation (above) transforms into the expression below.

<div align="center">

$H \Psi = E \Psi $

$\left[-\frac{1}{\gamma + 1} \nabla^2 + V\right] \Psi = E \Psi$

$\gamma = \sqrt{1 + \langle \Psi | \hat{p}^2 | \Psi \rangle / c^2}$

</div>

Below are the outputs of wavefunctions solved with this approximation.

> TODO: Paste graph of orbital contractions for each hydrogen-like ion across the periodic table. Compare to the first-order approximation on Wikipedia. Next, compare the calculated orbital energies to the exact solution from [this source](https://doi.org/10.1038/s41598-020-71505-w).
