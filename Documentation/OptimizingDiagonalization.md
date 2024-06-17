# Optimizating Diagonalization

This page discusses a proposed rewrite of the matrix diagonalizer, which would allow it to surpass Accelerate in performance.

## Tasks

Refactor the panel factorization into a separate file. That way,
it will be straightforward to optimize with BLAS and recursive panel
factorization. In addition, SBR will be easier to implement.
- Only accepts 'smallBlockSize' as a parameter. That way, it's the
  caller's responsibility to do the recursion part of recursive
  factorization, and to generate the T matrix in-place.
- Refactor the code to use T directly, without transposing it.

Employ the remaining improvements to single-core execution speed:
- exploiting symmetry
- eliding multiplications by 0
- consider successive band reductions (SBR)

Precompute the WY transforms and check that there are no major
regressions. Parallelize everything as much as possible:
- parallelize the last two stages
- parallelize the bulge chasing

The last few optimizations could be important in the future, but
perhaps not an economical use of time:
- pipeline the bulge chasing with panel factorization
- prepare the WY transforms simultaneously with other stages

## Benchmarking Code

Reporting diagnostics:

```swift
var bandReflectors: [Float] = []
var bulgeReflectors: [Float] = []

// Reduce the bandwidth from 'problemSize' to 'blockSize'.
let checkpoint0 = CACurrentMediaTime()
createMatrix(descriptor: descriptor)
bandReflectors = reduceToBandForm()

// If the matrix is already in tridiagonal form, there is no work to do.
let checkpoint1 = CACurrentMediaTime()
if blockSize > 1 {
  bulgeReflectors = chaseBulges()
}

// Acquire the eigenvectors (using the LAPACK divide-and-conquer until all
// other bottlenecks are suppressed).
let checkpoint2 = CACurrentMediaTime()
solveTridiagonalEigenproblem()

// If the matrix was already in tridiagonal form, there is no work to do.
let checkpoint3 = CACurrentMediaTime()
if blockSize > 1 {
  backTransform(bulgeReflectors: bulgeReflectors)
}

// Expand the bandwidth from 'blockSize' to 'problemSize'.
let checkpoint4 = CACurrentMediaTime()
backTransform(bandReflectors: bandReflectors)
destroyMatrix()
let checkpoint5 = CACurrentMediaTime()

// Report diagnostics for debugging performance.
if problemSize >= 100 {
  let time01 = 1e6 * (checkpoint1 - checkpoint0)
  let time12 = 1e6 * (checkpoint2 - checkpoint1)
  let time23 = 1e6 * (checkpoint3 - checkpoint2)
  let time34 = 1e6 * (checkpoint4 - checkpoint3)
  let time45 = 1e6 * (checkpoint5 - checkpoint4)
  let times = SIMD8(time01, time12, time23, time34, time45, 0, 0, 0)
  
  let total = times.sum()
  let percents = times / total * 100
  func format(_ number: Double) -> String {
    String(format: "%.1f", number)
  }
  func printPart(index: Int, label: String) {
    print("[\(format(percents[index]))%", terminator: ", ")
    print(format(times[index]), terminator: " μs] - ")
    print(label)
  }
  
  print()
  print("[n = \(problemSize)] Performance Breakdown:")
  printPart(index: 0, label: "Reduction to band form")
  printPart(index: 1, label: "Bulge chasing")
  printPart(index: 2, label: "Divide and conquer")
  printPart(index: 3, label: "Back transformation (2nd stage)")
  printPart(index: 4, label: "Back transformation (1st stage)")
}
```

Benchmarking successively larger matrix sizes:

```swift
// Compute-intensive benchmarks across a wide range of problem sizes.
benchmarkProblemSize(n: 125, trialCount: 3)
benchmarkProblemSize(n: 200, trialCount: 3)
benchmarkProblemSize(n: 250, trialCount: 3)
benchmarkProblemSize(n: 300, trialCount: 3)
benchmarkProblemSize(n: 400, trialCount: 3)
benchmarkProblemSize(n: 500, trialCount: 3)

// We need to speed up the custom eigensolver before testing these sizes.
benchmarkProblemSize(n: 600, trialCount: 3)
benchmarkProblemSize(n: 750, trialCount: 3)
benchmarkProblemSize(n: 800, trialCount: 3)
benchmarkProblemSize(n: 1000, trialCount: 3)
```

## Preliminary Benchmark Data

Implementations:
- Accelerate `dsyevd`
- Accelerate `ssyevd`
- Accelerate `ssyev`
- Current state of the `Diagonalize` API
  - With Swift release mode (`-Ounchecked`)
- Accelerate divide and conquer with the tridiagonal matrix from `Diagonalize`

Problems:
- Eigenvalues only, identity matrix (Λ1)
- Eigenvalues only, continuous spectrum (Λ2)
- Eigenvalues and eigenvectors, identity matrix (Σ1)
- Eigenvalues and eigenvectors, continuous spectrum (Σ2)

Benchmarking procedure:
- Problem sizes range from 147 to 407
- Per problem size: fastest of 3 samples
- Per 16-bin: fastest problem size

### Accelerate `dsyevd`

| n | Λ1 | Λ2 | Σ1 | Σ2 |
| - | -- | -- | -- | -- |
| 145&ndash;160 | 143.6 |   5.3 |   6.4 |   2.7 | 
| 161&ndash;176 | 167.8 |   6.2 |   6.9 |   3.1 | 
| 177&ndash;192 | 181.6 |   6.4 |   7.1 |   3.3 | 
| 193&ndash;208 | 198.0 |   6.5 |   7.3 |   3.6 | 
| 209&ndash;224 | 220.3 |   6.6 |   7.8 |   3.9 | 
| 225&ndash;240 | 228.0 |   6.7 |   7.8 |   4.0 | 
| 241&ndash;256 | 257.3 |   6.9 |   7.8 |   4.2 | 
| 257&ndash;272 | 272.0 |   6.9 |   7.8 |   4.4 | 
| 273&ndash;288 | 285.5 |   7.1 |   8.1 |   4.6 | 
| 289&ndash;304 | 309.7 |   7.3 |   8.3 |   4.8 | 
| 305&ndash;320 | 331.9 |   7.5 |   8.7 |   5.1 | 
| 321&ndash;336 | 348.4 |   7.5 |   8.8 |   5.1 | 
| 337&ndash;352 | 357.3 |   7.4 |   8.7 |   5.4 | 
| 353&ndash;368 | 379.5 |   7.6 |   8.7 |   5.4 | 
| 369&ndash;384 | 386.3 |   7.6 |   8.8 |   5.5 | 
| 385&ndash;400 | 405.6 |   8.0 |   9.5 |   5.9 | 
| 401&ndash;416 | 411.0 |   7.8 |   8.9 |   5.9 |

### Accelerate `ssyevd`

| n | Λ1 | Λ2 | Σ1 | Σ2 |
| - | -- | -- | -- | -- |
| 145&ndash;160 | 176.4 |   7.1 |   8.3 |   3.5 | 
| 161&ndash;176 | 198.9 |   7.7 |   8.5 |   3.8 | 
| 177&ndash;192 | 211.6 |   8.2 |   8.8 |   4.1 | 
| 193&ndash;208 | 225.6 |   8.6 |   9.1 |   4.5 | 
| 209&ndash;224 | 257.2 |   9.2 |  10.4 |   5.2 | 
| 225&ndash;240 | 276.3 |   9.5 |  10.4 |   5.6 | 
| 241&ndash;256 | 296.6 |   9.9 |  11.1 |   5.9 | 
| 257&ndash;272 | 316.2 |  10.0 |  11.3 |   6.1 | 
| 273&ndash;288 | 337.8 |  10.2 |  11.5 |   6.4 | 
| 289&ndash;304 | 354.7 |  10.5 |  11.9 |   6.8 | 
| 305&ndash;320 | 376.7 |  10.7 |  12.4 |   6.9 | 
| 321&ndash;336 | 397.1 |  10.9 |  12.4 |   7.1 | 
| 337&ndash;352 | 419.9 |  10.9 |  12.4 |   7.4 | 
| 353&ndash;368 | 437.4 |  11.1 |  12.9 |   8.0 | 
| 369&ndash;384 | 454.0 |  11.3 |  12.5 |   8.0 | 
| 385&ndash;400 | 474.0 |  11.9 |  13.6 |   8.6 | 
| 401&ndash;416 | 488.4 |  11.6 |  13.2 |   8.8 | 

### Accelerate `ssyev`

| n | Λ1 | Λ2 | Σ1 | Σ2 |
| - | -- | -- | -- | -- |
| 145&ndash;160 |  11.8 |   6.5 |  10.5 |   1.5 | 
| 161&ndash;176 |  12.5 |   7.2 |  10.2 |   1.6 | 
| 177&ndash;192 |  13.5 |   7.7 |  11.4 |   1.6 | 
| 193&ndash;208 |  14.0 |   8.2 |  11.2 |   1.7 | 
| 209&ndash;224 |  15.3 |   9.0 |  12.3 |   1.7 | 
| 225&ndash;240 |  15.4 |   9.3 |  12.4 |   1.7 | 
| 241&ndash;256 |  16.5 |   9.9 |  13.4 |   1.7 | 
| 257&ndash;272 |  16.2 |  10.2 |  12.8 |   1.7 | 
| 273&ndash;288 |  16.5 |  10.5 |  13.6 |   1.7 | 
| 289&ndash;304 |  17.0 |  10.9 |  13.6 |   1.7 | 
| 305&ndash;320 |  17.4 |  11.2 |  14.1 |   1.8 | 
| 321&ndash;336 |  17.6 |  11.5 |  14.2 |   1.8 | 
| 337&ndash;352 |  17.6 |  11.7 |  14.5 |   1.7 | 
| 353&ndash;368 |  17.8 |  12.0 |  14.6 |   1.7 | 
| 369&ndash;384 |  18.1 |  12.4 |  14.8 |   1.7 | 
| 385&ndash;400 |  19.3 |  13.1 |  15.8 |   1.7 | 
| 401&ndash;416 |  18.4 |  12.8 |  14.9 |   1.7 | 

### Current state of `Diagonalize`

| n | Σ1 | Σ2 |
| - | -- | -- |
| 145&ndash;160 |   2.8 |   2.1 | 
| 161&ndash;176 |   3.0 |   2.2 | 
| 177&ndash;192 |   3.3 |   2.4 | 
| 193&ndash;208 |   3.3 |   2.6 | 
| 209&ndash;224 |   3.9 |   2.9 | 
| 225&ndash;240 |   3.8 |   3.0 | 
| 241&ndash;256 |   4.4 |   3.4 | 
| 257&ndash;272 |   4.6 |   3.5 | 
| 273&ndash;288 |   4.9 |   3.8 | 
| 289&ndash;304 |   5.0 |   3.8 | 
| 305&ndash;320 |   5.2 |   4.1 | 
| 321&ndash;336 |   5.3 |   4.2 | 
| 337&ndash;352 |   5.1 |   4.3 | 
| 353&ndash;368 |   5.4 |   4.5 | 
| 369&ndash;384 |   5.7 |   4.5 | 
| 385&ndash;400 |   6.2 |   4.9 | 
| 401&ndash;416 |   6.2 |   4.9 | 

### Accelerate divide and conquer

| n | Λ1 | Λ2 |
| - | -- | -- |
| 145&ndash;160 | 148.9 |   6.9 | 
| 161&ndash;176 | 185.0 |   8.2 | 
| 177&ndash;192 | 198.7 |   8.9 | 
| 193&ndash;208 | 220.7 |  10.1 | 
| 209&ndash;224 | 238.8 |  11.3 | 
| 225&ndash;240 | 256.1 |  12.3 | 
| 241&ndash;256 | 274.4 |  13.6 | 
| 257&ndash;272 | 292.4 |  14.6 | 
| 273&ndash;288 | 307.0 |  15.6 | 
| 289&ndash;304 | 328.2 |  16.9 | 
| 305&ndash;320 | 342.6 |  18.0 | 
| 321&ndash;336 | 353.6 |  19.5 | 
| 337&ndash;352 | 354.6 |  20.5 | 
| 353&ndash;368 | 371.6 |  21.9 | 
| 369&ndash;384 | 383.1 |  23.1 | 
| 385&ndash;400 | 401.3 |  24.5 | 
| 401&ndash;416 | 434.0 |  25.4 | 

| n | Λ1 | Λ2 |
| - | -- | -- |
| 369&ndash;384 | 377.6 |  23.4 | 
| 433&ndash;448 | 465.9 |  28.9 | 
| 497&ndash;512 | 472.7 |  34.5 | 
| 561&ndash;576 | 588.1 |  41.0 | 
| 625&ndash;640 | 534.9 |  46.7 | 
| 689&ndash;704 | 587.4 |  53.0 | 
| 753&ndash;768 | 636.2 |  59.9 | 
| 817&ndash;832 | 585.3 |  62.7 | 
| 881&ndash;896 | 737.8 |  73.3 | 

With a 407x407 hamiltonian, we can expect the following latencies for 7 SCF cycles:
- Accelerate `dsyevd`: 80 ms
- Accelerate `ssyevd`: 54 ms
- If divide and conquer is the bottleneck: 19 ms

With an 896x896 hamiltonian, we can expect the following latencies for 7 SCF cycles:
- Accelerate `dsyevd`: 629 ms
- Accelerate `ssyevd`: 379 ms
- If divide and conquer is the bottleneck: 69 ms

We should prove that the bottleneck in xTB is actually matrix diagonalization at these problem sizes. What percent of the execution time would be spent on diagonalization, with only the latency of divide and conquer?

| Problem Size | SCF Cycles | SCF Latency (Reported) | Total Latency (Reported) | Actual Latency |
| ------------ | ---------- | ---------------------- | ------------------------ | --- |
| 1 | 2 | 3 | 4 | 5 |
