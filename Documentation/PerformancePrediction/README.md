# Performance Prediction

This experiment will predict the performance of matrix diagonalization with different algorithms and hardware. There is a Swift script, which takes different parameters like block size, computer power, parallelism, etc. If written correctly, it should reproduce the performance behavior of LAPACK with one-stage tridiagonalization.

> The divide-and-conquer part will be omitted from the performance prediction. It is out of scope, and not the bottleneck being optimized.

## QR Factorization

Predictions of execution speed of recursively blocked QR factorization, with various matrix sizes. For rectangular matrices, GFLOPS/k is defined with $O(k \cdot nb(ns)^2)$. $nb$ is the dimension of the large block, and $ns$ is the dimension of the small block. This notation is relevant because SBR switches between block sizes during each stage.

I'll compare the slowness of QR factorization to the cost of the trailing matrix updates. If it dominates execution time for target problem sizes, I might need to analyze TSQR. TSQR employs many small QR factorizations, which would be simulated well with this section's experiment.

> NOTE: TSQR could utilize multiple CPU cores, if the overall latency exceeds 10 microseconds.

Below, I simulate various hierarchical blocking schemes on square matrices.

### 1-8-32 Blocking Scheme

| n    | Time (1x1) | Time (8x8) | Time (32x32) | Latency | GFLOPS/k |
| ---- | ---------- | ---------- | ------------ | ------- | -------- |
| 64   | 51.1%      | 43.6%      | 5.3%         | 8 μs    | 35       |
| 128  | 48.6%      | 42.1%      | 9.4%         | 30 μs   | 69       |
| 256  | 44.6%      | 38.9%      | 16.5%        | 129 μs  | 130      |
| 512  | 38.5%      | 33.7%      | 27.8%        | 589 μs  | 228      |
| 1024 | 30.3%      | 26.6%      | 43.2%        | 2.98 ms | 360      |
| 2048 | 21.2%      | 18.7%      | 60.1%        | 17.0 ms | 507      |
| asymptote | 0.0%  | 0.0%       | 100.0%       | -       | 850      |

<details>
<summary>Code</summary>

```swift
let argumentN = CommandLine.arguments[1]
let n: Float = Float(argumentN)!

let latencyForHouseholder1x1: Float = 2 * n * (1 * 1) / 12e9

let latencyForT8x8: Float = n * (8 * 8) / 42.5e9 + (8 * 8 * 8) / 12e9

let latencyForPanel8x8: Float = 0.5 * (8 * 8) * latencyForHouseholder1x1 + latencyForT8x8

let latencyForVTV8x8: Float = 2 * n * (8 * 8) / 42.5e9 + (8 * 8 * 8) / 42.5e9

let latencyForPanel32x32: Float = (32 / 8) * latencyForPanel8x8 + 0.5 * (32 / 8) * (32 / 8) * latencyForVTV8x8

// Assume the AMX is writing out of registers, with a cache bandwidth of 400 GB/s.
let latencyToWrite32x32Matrix: Float = (32 * 32 * 4) / 400e9

// Assume the 32x32 panel multiplies by a block-diagonal T matrix.
let latencyForVTV32x32: Float = 2 * n * (32 * 32) / 566.2e9 + (32 * 32 * 32) / 566.2e9 + latencyToWrite32x32Matrix

print("latency to apply individual reflector:", latencyForHouseholder1x1)
print("latency to generate T:", latencyForT8x8)
print("latency to construct 8x8 panel:", latencyForPanel8x8)
print("latency to apply 8x8 reflector block:", latencyForVTV8x8)
print("latency to construct 32x32 panel:", latencyForPanel32x32)
print("latency to write 32x32 matrix:", latencyToWrite32x32Matrix)
print("latency to apply 32x32 reflector block:", latencyForVTV32x32)

print("GFLOPS/k to apply individual reflector:", n * (1 * 1) / latencyForHouseholder1x1 / 1e9)
print("GFLOPS/k to generate T:", (n * 8 * 8) / latencyForT8x8 / 1e9)
print("GFLOPS/k to construct 8x8 panel:", n * (8 * 8) / latencyForPanel8x8 / 1e9)
print("GFLOPS/k to apply 8x8 reflector block:", n * (8 * 8) / latencyForVTV8x8 / 1e9)
print("GFLOPS/k to construct 32x32 panel:", n * (32 * 32) / latencyForPanel32x32 / 1e9)
print("GFLOPS/k to apply 32x32 reflector block:", n * (32 * 32) / latencyForVTV32x32 / 1e9)

let timeForPanel8x8: Float = (n / 8) * latencyForPanel8x8
let timeForVTV8x8: Float = 0.5 * (n / 8) * (32 / 8) * latencyForVTV8x8
let timeForVTV32x32: Float = 0.333 * (n / 32) * (n / 32) * latencyForVTV32x32
let timeTotal = timeForPanel8x8 + timeForVTV8x8 + timeForVTV32x32

print("time spent constructing panels (8x8):", 100 * timeForPanel8x8 / timeTotal)
print("time spent applying panels (8x8):", 100 * timeForVTV8x8 / timeTotal)
print("time spent applying panels (32x32):", 100 * timeForVTV32x32 / timeTotal)
print("latency for entire matrix:", timeTotal * 1e6, "μs")
print("GFLOPS/k for entire matrix:", n * n * n / timeTotal / 1e9)
print("GFLOPS/k for same size GEMM:", 566.2)
```

</details>

<details>
<summary>Output for n=1024</summary>

```
latency to apply individual reflector: 1.7066667e-07
latency to generate T: 1.5846902e-06
latency to construct 8x8 panel: 7.0460237e-06
latency to apply 8x8 reflector block: 3.0960941e-06
latency to construct 32x32 panel: 5.2952848e-05
latency to write 32x32 matrix: 1.024e-08
latency to apply 32x32 reflector block: 3.7720206e-06
GFLOPS/k to apply individual reflector: 6.0
GFLOPS/k to generate T: 41.35572
GFLOPS/k to construct 8x8 panel: 9.301132
GFLOPS/k to apply 8x8 reflector block: 21.167315
GFLOPS/k to construct 32x32 panel: 19.802069
GFLOPS/k to apply 32x32 reflector block: 277.98788
time spent constructing panels (8x8): 30.25749
time spent applying panels (8x8): 26.590893
time spent applying panels (32x32): 43.15162
latency for entire matrix: 2980.72 μs
GFLOPS/k for entire matrix: 360.229
GFLOPS/k for same size GEMM: 566.2
```

</details>

### Latency to Factorize an 8x8 Tile

Here, I predict the latency to factorize an 8x8 tile on different processors. The analysis assumes a 1-8 blocking scheme.

From handwritten calculations, I should expect the following for `n = 1024`. `1x1` is the latency to generate a single reflector, while `8x8` is the latency to factorize the entire panel. The 1/2 constant factor preceding the application of 1x1 reflectors is omitted, to keep the calculations consistent with each other.

| Processor                 | Latency (1x1) | Latency (8x8, LL) | Latency (8x8, RL) |
| ------------------------- | ------------- | ----------------- | ----------------- |
| CPU                       | 72 ns         | 7.2 μs            | 6.1 μs            |
| CPU (AMX accelerated)     | 72 ns         | 5.3 μs            | 4.2 μs            |
| GPU (serial over SIMDs)   | 228 ns        | 5.9 μs            | 5.9 μs            |
| GPU (parallel over SIMDs) | 135 ns        | 8.3 μs            | 3.5 μs            |

For `n = 1024`, the latencies for GPU and AMX theoretically outperform that of a CPU core. However, both are bound by a large `O(1)` latency. The GPU must synchronize over all SIMDs in a threadgroup, which is assumed to take [86 ns](https://chipsandcheese.com/2022/05/21/igpu-cache-setups-compared-including-m1/). The CPU must incur latency to switch contexts to the AMX. The minimum overhead for GPU (parallel over SIMDs) is 1.38 μs.

Next, the above calculations are reproduced in code, allowing them to be re-run for different values of `n`.

| n = 256                   | Latency (1x1) | Latency (8x8, LL) | Latency (8x8, RL) |
| ------------------------- | ------------- | ----------------- | ----------------- |
| CPU                       |
| CPU (AMX accelerated)     |
| GPU (serial over SIMDs)   | 111 ns        | 2.0 μs            | 2.0 μs            |
| GPU (parallel over SIMDs) |

| n = 512                   | Latency (1x1) | Latency (8x8, LL) | Latency (8x8, RL) |
| ------------------------- | ------------- | ----------------- | ----------------- |
| CPU                       |
| CPU (AMX accelerated)     |
| GPU (serial over SIMDs)   |
| GPU (parallel over SIMDs) | 148 ns        | 3.3 μs            | 3.3 μs            |

| n = 1024                  | Latency (1x1) | Latency (8x8, LL) | Latency (8x8, RL) |
| ------------------------- | ------------- | ----------------- | ----------------- |
| CPU                       |
| CPU (AMX accelerated)     |
| GPU (serial over SIMDs)   | 228 ns        | 5.9 μs            | 5.9 μs            |
| GPU (parallel over SIMDs) |

| n = 2048                  | Latency (1x1) | Latency (8x8, LL) | Latency (8x8, RL) |
| ------------------------- | ------------- | ----------------- | ----------------- |
| CPU                       |
| CPU (AMX accelerated)     |
| GPU (serial over SIMDs)   | 426 ns        | 11.4 μs           | 11.4 μs           |
| GPU (parallel over SIMDs) |

### 1-2-8-32 Blocking Scheme
