# Relativity Corrections

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
