# Dense Diagonalization Experiment

In this experiment, I attempt to create an algorithm that diagonalizes a dense
matrix. Every part of the algorithm can be perfectly parallelized.

```
// Total cost: 1 matrix multiplication
parallel-power-iteration(H, inout Ψ)
  HΨ = H * Ψ       O(n^3)
  Ψ <- HΨ / ||HΨ|| O(n^2)

// Total cost: 1 matrix multiplication
parallel-steepest-descent-iteration(H, inout Ψ, inout HΨ)
  E = vector-wise Ψ * HΨ  O(n^2)
  rayleigh = E / ||Ψ||    O(n^2)
  r = HΨ - rayleigh * Ψ   O(n^2)
  yield max ||r||         O(n^2)
  
  Hr = H * r              O(n^3)
  construct mm (          O(n^2)
    Ψ, HΨ, r, Hr, E)
  construct λ (m)         O(n)
  if λ[i] == 0
    lock vector i

  Ψ <- Ψ + λr             O(n^2)
  HΨ <- HΨ + λHr          O(n^2)

// Total cost: 3 symmetric matrix multiplications
parallel-orthogonalize-iteration(inout Ψ)
  W = Ψ^T * Ψ             O(n^3)
  yield max |W|           O(n^2)

  W <- lower triangle (W) O(n^2)
  F = -W * Ψ              O(n^3)
  l = vector-wise ||F||   O(n^2)
  construct λ (l)         O(n^2)
  
  F = -λW * Ψ             O(n^3)
  Ψ <- Ψ + F              O(n^2)
```

Here are some high-level solvers constructed from the above primitives.

```
inertia-tensor-solver(H)
  E = factor-characteristic-polynomial(H) O(n^3)
  Ψ = identity matrix                     O(n^2)
  repeat until ΨHΨ == E {
    parallel-power-iteration(H, Ψ)        O(n^3)
    parallel-orthogonalize-iteration(Ψ)   O(n^3)
  }
  return (E, Ψ)

self-consistent-field-solver()
  Ψ = ansatz                               O(n^2)
  repeat until ||r|| < 1e-7 {
    H = construct hamiltonian(Ψ)           O(n^2)
    HΨ = H * Ψ                             O(n^3)
    repeat 5 times
      parallel-steepest-descent-iteration( O(n^3)
        H, Ψ, HΨ)
    
    Ψ = normalize(Ψ)                       O(n^2)
    repeat until ||W|| < 1e-7
      parallel-orthogonalize-iteration(Ψ)  O(n^3)
      Ψ = normalize(Ψ)                     O(n^2)
  }
  return Ψ

simultaneous-dense-solver(H, ν1, ν2)
  Ψ = identity matrix                      O(n^2)
  repeat until ||r|| < 1e-7,
               ||W|| < 1e-7 {
    HΨ = H * Ψ                             O(n^3)
    repeat ν1 times
      parallel-steepest-descent-iteration( O(n^3)
        H, Ψ, HΨ)
    Ψ = normalize(Ψ)                       O(n^2)
    
    repeat ν2 times
      parallel-orthogonalize-iteration(Ψ)  O(n^3)
      Ψ = normalize(Ψ)                     O(n^2)
  }
  E = Ψ * H * Ψ                            O(n^3)
  return (E, Ψ)
```

Next, we will test the solvers on some artificially generated, diagonalizable
matrices. Input a set of manually selected eigenvalues, randomly generate a
set of orthonormal eigenvectors, then construct the test matrix. The ultimate
goal is to benchmark GFLOPS/k, using the Apple AMX for matrix multiplications.
Compare the custom algorithm to the LAPACK diagonalizer.

```
GFLOPS/k = (matrix dimension)^3 / (time to solution) / 1e9
```

## Orthonormalization

I found that the "fast orthogonalization" algorithm above had some flaws. There seems to be a critical slowing down, where each iteration only annihiliates one column. I settled for a panel-factorized version of Gram-Schmidt. This should have near-optimal performance, and is still amenable to parallelization within a single compute node.

I also created a custom Metal compiler that supports async copies. There will be a need to author dense linear algebra kernels for the Apple GPU.
