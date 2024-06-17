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
