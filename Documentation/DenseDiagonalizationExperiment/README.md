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

I also created a custom Metal compiler that supports async copies. There may be a need to author dense linear algebra kernels for the Apple GPU.

## Diagonalization

The attempts at a custom iterative diagonalizer failed. Therefore, I am relying on the LAPACK implementation for now.

Here is the `GFLOPS/k` for various methods of decomposition. Two-stage tridiagonalization fixes the bottleneck, but the results are often inaccurate. Using double precision does not fix the results. That suggests it could be a bug in Accelerate or the LAPACK open-source reference code.

| n   | Algorithm | Λ | Σ | (Λ, Σ) |
| --: | :-------- | --: | --: | --: |
| 1000 | `ssyev_`         | 18.7 | 3.0 | 2.6 |
| 1000 | `ssyevd_`        | 14.8 | 130.4 | 13.3 |
| 1000 | `ssyevd_2stage_` | 60.5 | -   | -   |

| n   | Algorithm | Λ | Σ | (Λ, Σ) |
| --: | :-------- | --: | --: | --: |
|  125 | `ssyev_`         | 5.2 | 0.8 | 0.7 |
|  250 | `ssyev_`         | 9.8 | 1.9 | 1.6 |
|  500 | `ssyev_`         | 15.0 | 1.9 | 1.7 |
|  750 | `ssyev_`         | 19.2 | 2.6 | 2.3 |
| 1000 | `ssyev_`         | 18.7 | 3.0 | 2.6 |
|  125 | `ssyevd_`        | 5.8 | 5.3 | 2.8 |
|  250 | `ssyevd_`        | 9.7 | 14.2 | 5.8 |
|  500 | `ssyevd_`        | 13.3 | 42.1 | 10.1 |
|  750 | `ssyevd_`        | 15.6 | 66.9 | 12.6 |
| 1000 | `ssyevd_`        | 14.8 | 130.4 | 13.3 |
|  125 | `ssyevd_2stage_` | 5.0 | -   | -   |
|  250 | `ssyevd_2stage_` | 11.2 | -   | -   |
|  500 | `ssyevd_2stage_` | 26.3 | -   | -   |
|  750 | `ssyevd_2stage_` | 38.0 | -   | -   |
| 1000 | `ssyevd_2stage_` | 60.5 | -   | -   |
| 1250 | `ssyevd_2stage_` | 77.4 | -   | -   |
| 1500 | `ssyevd_2stage_` | 79.6 | -   | -   |

> Table of GFLOPS/k performance. The maximum possible is ~750-1500 for Apple AMX with GEMM. Most data was truncated at n ≤ 1000 because benchmark quality deteriorated rapidly for larger matrices.

It seems worthwhile to write a two-stage kernel for the AMX from scratch. That should both fix the bug, and allow `ssyevd_2stage_` to generate eigenvectors.

## Accelerate Performance

I am deciding whether to attempt to harness the AMX instructions from Swift, or use the built-in GEMM kernel from Accelerate.

> The tables show maximum GFLOPS/k on each processor.

| Block Size | NN     | NT    | TN    | TT    | Best   |
| :----------: | :------: | :-----: | :-----: | :-----: | :-----: |
| 1x1        | 1.6    | 1.6   | 1.6   | 1.6   | 1.6    |
| 2x2        | 3.8    | 1.6   | 5.7   | 1.8   | 5.7    |
| 4x4        | 12.4   | 12.3  | 18.2  | 12.3  | 18.2   |
| 8x8        | 27.7   | 42.5  | 23.1  | 28.7  | 42.5   |
| 16x16      | 98.3   | 172.4 | 71.5  | 95.3  | 172.4  |
| 32x32      | 267.8  | 566.2 | 172.4 | 262.2 | 566.2  |
| 64x64      | 377.9  | 581.9 | 290.1 | 378.5 | 581.9  |
| 128x128    | 499.7  | 592.0 | 408.3 | 499.9 | 592.0  |
| 256x256    | 995.4  | 989.9 | 733.9 | 831.7 | 995.4  |
| 512x512    | 1071.1 | 892.8 | 876.3 | 952.5 | 1071.1 |
| 1024x1024  | 917.3  | 772.5 | 942.3 | 885.3 | 942.3  |

| Block Size | Best   | AMX Theoretical Limit | AMX Utilization |
| :----------: | :------: | :---------------------: | :---------------: |
| 1x1        | 1.6    | 0.7                   | 228.6%          |
| 2x2        | 5.7    | 2.9                   | 196.6%          |
| 4x4        | 18.2   | 12                    | 151.7%          |
| 8x8        | 42.5   | 46                    | 92.4%           |
| 16x16      | 172.4  | 183                   | 94.2%           |
| 32x32      | 566.2  | 733                   | 77.2%           |
| 64x64      | 581.9  | 740                   | 78.6%           |
| 128x128    | 592.0  | 1470                  | 40.3%           |
| 256x256    | 995.4  | 1480                  | 67.3%           |
| 512x512    | 1071.1 | 1480                  | 72.4%           |
| 1024x1024  | 942.3  | 1480                  | 63.7%           |

| Block Size | NEON | Accelerate | GPU       | Fastest    |
| :----------: | :----: | :----------: | :---------: | :----------: |
| 1x1        | 25   | 1.6        | 16        | NEON       |
| 2x2        | 50   | 5.7        | 4–21      | NEON       |
| 4x4        | 50   | 18.2       | 8–42      | NEON       |
| 8x8        | 50   | 42.5       | 17–83     | Accelerate |
| 16x16      | 50   | 172.4      | 66–83     | Accelerate |
| 32x32      | 105  | 566.2      | 175–262   | Accelerate |
| 64x64      | 400  | 581.9      | 175–524   | Accelerate |
| 128x128    | 400  | 592.0      | 466–1400  | Accelerate |
| 256x256    | 400  | 995.4      | 1500–2560 | GPU        |
| 512x512    | 400  | 1071.1     | 3000      | GPU        |
| 1024x1024  | 400  | 942.3      | 3800      | GPU        |

Conclusion:
- NEON (Swift SIMD) is suitable for block sizes 1&mdash;8.
- AMX (Accelerate `sgemm_`) is suitable for block sizes 8&mdash;128.
- GPU (Metal `simdgroup_async_copy`) is suitable for block sizes 128&mdash;1024.
