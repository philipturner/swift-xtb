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

<details>
<summary>Code</summary>

```swift
// MARK: - Functions

func gpuSerial(n: Int) {
  print()
  print("GPU serial")
  
  // Generating the reflector.
  var cyclesGenerate: Int
  do {
    let cyclesZeroOut: Int = max((n / 32) * 2, 4)
    let cyclesFMADot: Int = max((n / 32) * 2, 4)
    let cyclesShuffle: Int = 5 * 4
    let cyclesRsqrt: Int = 8
    let cyclesDiv: Int = 2 * 6
    let cyclesFMAScale: Int = max((n / 32) * 2, 4)
    let cyclesWrite: Int = max((n / 32) * 2, 56)
    print("zero out  ", cyclesZeroOut)
    print("FMA dot   ", cyclesFMADot)
    print("shuffle   ", cyclesShuffle)
    print("RSQRT     ", cyclesRsqrt)
    print("serial DIV", cyclesDiv)
    print("FMA scale ", cyclesFMAScale)
    print("write     ", cyclesWrite)
    
    cyclesGenerate =
    cyclesZeroOut + cyclesFMADot + cyclesShuffle + cyclesRsqrt +
    cyclesDiv + cyclesFMAScale + cyclesWrite
    print("generate reflector (total)", cyclesGenerate)
  }
  
  // Applying the reflector.
  var cyclesApply: Int
  do {
    let cyclesLoad: Int = max(8 * n / 16, 56)
    let cyclesFMADot: Int = max((n / 32) * 2, 4)
    let cyclesShuffle: Int = 5 * 4
    let cyclesFMAScale: Int = max((n / 32) * 2, 4)
    print("load (core)      ", cyclesLoad)
    print("FMA (core)       ", cyclesFMADot)
    print("shuffle (latency)", cyclesShuffle)
    print("FMA (core)       ", cyclesFMAScale)
    
    cyclesApply =
    cyclesLoad + cyclesFMADot + cyclesShuffle + cyclesFMAScale
    print("apply reflector (total)", cyclesApply)
  }
  
  // Repeat the above subroutines 8 times.
  let cycles8x8 = 8 * (cyclesGenerate + cyclesApply)
  let latency1x1 = Float(cyclesGenerate) / 1.296e9
  let latency8x8 = Float(cycles8x8) / 1.296e9
  print("latency (1x1)", latency1x1)
  print("latency (8x8)", latency8x8)
}

func gpuParallel(n: Int, rightLooking: Bool) {
  print()
  print("GPU parallel", terminator: " ")
  if rightLooking {
    print("(RL)")
  } else {
    print("(LL)")
  }
  
  // Generating the reflector.
  var cyclesGenerate: Int
  do {
    let cyclesZeroOut: Int = max(n / 256 * 2, 4)
    let cyclesFMADot: Int = max(n / 256 * 2, 4)
    let cyclesShuffle: Int = 5 * 4
    let cyclesBroadcast: Int = 2 * 56
    let cyclesRsqrt: Int = 8
    let cyclesDiv: Int = 2 * 6
    let cyclesFMAScale: Int = max(n / 256 * 2, 4)
    print("zero out  ", cyclesZeroOut)
    print("FMA dot   ", cyclesFMADot)
    print("shuffle   ", cyclesShuffle)
    print("broadcast ", cyclesBroadcast)
    print("RSQRT     ", cyclesRsqrt)
    print("serial DIV", cyclesDiv)
    print("FMA scale ", cyclesFMAScale)
    
    cyclesGenerate =
    cyclesZeroOut + cyclesFMADot + cyclesShuffle + cyclesBroadcast +
    cyclesRsqrt + cyclesDiv + cyclesFMAScale
    print("generate reflector (total)", cyclesGenerate)
  }
  
  // Applying the reflector.
  var cyclesApply: Int
  if rightLooking {
    let cyclesFMADot: Int = 8 * max(n / 256 * 2, 4)
    let cyclesShuffle: Int = 8 * 5 * 4
    let cyclesBroadcast: Int = 2 * 56
    let cyclesFMAScale: Int = 8 * max(n / 256 * 2, 4)
    print("FMA      ", cyclesFMADot)
    print("shuffle  ", cyclesShuffle)
    print("broadcast", cyclesBroadcast)
    print("FMA      ", cyclesFMAScale)
    
    cyclesApply =
    cyclesFMADot + cyclesShuffle + cyclesBroadcast + cyclesFMAScale
    print("apply reflector (total)", cyclesApply)
  } else {
    let cyclesFMADot: Int = max(n / 256 * 2, 4)
    let cyclesShuffle: Int = 5 * 4
    let cyclesBroadcast: Int = 2 * 56
    let cyclesFMAScale: Int = max(n / 256 * 2, 4)
    print("FMA      ", cyclesFMADot)
    print("shuffle  ", cyclesShuffle)
    print("broadcast", cyclesBroadcast)
    print("FMA      ", cyclesFMAScale)
    
    let cycles1 =
    cyclesFMADot + cyclesShuffle + cyclesBroadcast + cyclesFMAScale
    print("1 reflector", cycles1)
    
    cyclesApply = 8 * cycles1
    print("apply reflector (total)", cyclesApply)
  }
  
  // Repeat the above subroutines 8 times.
  let cycles8x8 = 8 * (cyclesGenerate + cyclesApply)
  let latency1x1 = Float(cyclesGenerate) / 1.296e9
  let latency8x8 = Float(cycles8x8) / 1.296e9
  print("latency (1x1)", latency1x1)
  print("latency (8x8)", latency8x8)
}

func cpu(n: Int, useAMX: Bool, rightLooking: Bool) {
  print()
  print("CPU", terminator: " ")
  if useAMX {
    print("AMX", terminator: " ")
  }
  if rightLooking {
    print("(RL)")
  } else {
    print("(LL)")
  }
  
  // Generating the reflector.
  var cyclesGenerate: Int
  do {
    let cyclesReadDot: Int = max(n / 16, 4 + 3)
    let cyclesReduce: Int = min(n / 16, 64).trailingZeroBitCount * 3
    let cyclesRsqrt: Int = 10
    let cyclesDiv: Int = 2 * 8
    let cyclesReadScale: Int = max(n / 16, 4 + 3)
    let cyclesWrite: Int = max(n / 16, 4)
    //    print("read dot  ", cyclesReadDot)
    //    print("reduce    ", cyclesReduce)
    //    print("RSQRT     ", cyclesRsqrt)
    //    print("serial DIV", cyclesDiv)
    //    print("read scale", cyclesReadScale)
    //    print("write     ", cyclesWrite)
    
    cyclesGenerate =
    cyclesReadDot + cyclesReduce + cyclesRsqrt + cyclesDiv +
    cyclesReadScale + cyclesWrite
    //    print("generate reflector (total)", cyclesGenerate)
  }
  
  // Applying the reflector.
  var cyclesApply: Int
  if useAMX {
    if rightLooking {
      let cyclesReadV1: Int = max(n / 32, 4)
      let cyclesReadA1: Int = 8 * max(n / 32, 4)
      let cyclesFMADot: Int = 8 * max(n / 32, 4)
      let cyclesReduce: Int = 8 * min(n / 32, 64).trailingZeroBitCount * 4
      let cyclesReadV2: Int = max(n / 32, 4)
      let cyclesReadA2: Int = 8 * max(n / 32, 4)
      let cyclesFMAScale: Int = 8 * max(n / 32, 4)
      let cyclesWrite: Int = 8 * max(n / 32, 4)
      //      print("read V (1)", cyclesReadV1)
      //      print("read A (1)", cyclesReadA1)
      //      print("FMA dot   ", cyclesFMADot)
      //      print("reduce    ", cyclesReduce)
      //      print("read V (2)", cyclesReadV2)
      //      print("read A (2)", cyclesReadA2)
      //      print("FMA scale ", cyclesFMAScale)
      //      print("write     ", cyclesWrite)
      
      cyclesApply =
      cyclesReadV1 + cyclesReadA1 + cyclesFMADot + cyclesReduce +
      cyclesReadV2 + cyclesReadA2 + cyclesFMAScale + cyclesWrite
      //      print("apply reflector (total)", cyclesApply)
    } else {
      let cyclesReadV1: Int = max(n / 32, 4)
      let cyclesReadA1: Int = max(n / 32, 4)
      let cyclesFMADot: Int = max(n / 32, 4)
      let cyclesReduce: Int = min(n / 32, 64).trailingZeroBitCount * 4
      let cyclesReadV2: Int = max(n / 32, 4)
      let cyclesReadA2: Int = max(n / 32, 4)
      let cyclesFMAScale: Int = max(n / 32, 4)
      let cyclesWrite: Int = max(n / 32, 4)
      //      print("read V (1)", cyclesReadV1)
      //      print("read A (1)", cyclesReadA1)
      //      print("FMA dot   ", cyclesFMADot)
      //      print("reduce    ", cyclesReduce)
      //      print("read V (2)", cyclesReadV2)
      //      print("read A (2)", cyclesReadA2)
      //      print("FMA scale ", cyclesFMAScale)
      //      print("write     ", cyclesWrite)
      
      let cycles1 =
      cyclesReadV1 + cyclesReadA1 + cyclesFMADot + cyclesReduce +
      cyclesReadV2 + cyclesReadA2 + cyclesFMAScale + cyclesWrite
      //      print("1 reflector", cycles1)
      
      cyclesApply = 8 * cycles1
      //      print("apply reflector (total)", cyclesApply)
    }
  } else {
    if rightLooking {
      let cyclesReadV1: Int = max(n / 16, 4)
      let cyclesReadA1: Int = 8 * max(n / 16, 4 + 3)
      let cyclesReduce: Int = 8 * min(n / 16, 64).trailingZeroBitCount * 3
      let cyclesReadV2: Int = max(n / 16, 4)
      let cyclesReadA2: Int = 8 * max(n / 16, 4 + 3)
      let cyclesWrite: Int = 8 * max(n / 16, 4)
      //      print("read V (1)", cyclesReadV1)
      //      print("read A (1)", cyclesReadA1)
      //      print("reduce    ", cyclesReduce)
      //      print("read V (2)", cyclesReadV2)
      //      print("read A (2)", cyclesReadA2)
      //      print("write     ", cyclesWrite)
      
      cyclesApply =
      cyclesReadV1 + cyclesReadA1 + cyclesReduce + cyclesReadV2 +
      cyclesReadA2 + cyclesWrite
      //      print("apply reflector (total)", cyclesApply)
    } else {
      let cyclesReadV1: Int = max(n / 16, 4)
      let cyclesReadA1: Int = max(n / 16, 4 + 3)
      let cyclesReduce: Int = min(n / 16, 64).trailingZeroBitCount * 3
      let cyclesReadV2: Int = max(n / 16, 4)
      let cyclesReadA2: Int = max(n / 16, 4 + 3)
      let cyclesWrite: Int = max(n / 16, 4)
      //      print("read V (1)", cyclesReadV1)
      //      print("read A (1)", cyclesReadA1)
      //      print("reduce    ", cyclesReduce)
      //      print("read V (2)", cyclesReadV2)
      //      print("read A (2)", cyclesReadA2)
      //      print("write     ", cyclesWrite)
      
      let cycles1 =
      cyclesReadV1 + cyclesReadA1 + cyclesReduce + cyclesReadV2 +
      cyclesReadA2 + cyclesWrite
      //      print("1 reflector", cycles1)
      
      cyclesApply = 8 * cycles1
      //      print("apply reflector (total)", cyclesApply)
    }
  }
  
  
  // Repeat the above subroutines 8 times.
  let cycles8x8 = 8 * (cyclesGenerate + cyclesApply)
  let latency1x1 = Float(cyclesGenerate) / 3.228e9
  let latency8x8 = Float(cycles8x8) / 3.228e9
  print("latency (1x1)", latency1x1)
  print("latency (8x8)", latency8x8)
}

// MARK: - Script

let argumentN = CommandLine.arguments[1]
let n: Int = Int(argumentN)!
cpu(n: n, useAMX: false, rightLooking: false)
cpu(n: n, useAMX: false, rightLooking: true)
cpu(n: n, useAMX: true, rightLooking: false)
cpu(n: n, useAMX: true, rightLooking: true)
```

</details>

| CPU                       | Latency (1x1) | Latency (8x8, LL) | Latency (8x8, RL) |
| ------------------------- | ------------- | ----------------- | ----------------- |
| n = 32                    | 15 ns         | 690 ns            | 550 ns            |
| n = 64                    | 15 ns         | 750 ns            | 620 ns            |
| n = 128                   | 18 ns         | 1.1 μs            | 840 ns            |
| n = 256                   | 27 ns         | 2.0 μs            | 1.5 μs            |
| n = 512                   | 42 ns         | 3.8 μs            | 2.7 μs            |
| n = 1024                  | 73 ns         | 7.2 μs            | 5.1 μs            |
| n = 2048                  | 133 ns        | 14.1 μs           | 9.7 μs            |

| CPU (AMX accelerated)     | Latency (1x1) | Latency (8x8, LL) | Latency (8x8, RL) |
| ------------------------- | ------------- | ----------------- | ----------------- |
| n = 32                    | 15 ns         | 670 ns            | 530 ns            |
| n = 64                    | 15 ns         | 760 ns            | 620 ns            |
| n = 128                   | 18 ns         | 860 ns            | 720 ns            |
| n = 256                   | 27 ns         | 1.6 μs            | 1.3 μs            |
| n = 512                   | 42 ns         | 2.9 μs            | 2.3 μs            |
| n = 1024                  | 73 ns         | 5.4 μs            | 4.3 μs            |
| n = 2048                  | 133 ns        | 10.4 μs           | 8.2 μs            |

| GPU (serial over SIMDs)   | Latency (1x1) | Latency (8x8, LL) | Latency (8x8, RL) |
| ------------------------- | ------------- | ----------------- | ----------------- |
| n = 32                    | 83 ns         | 1.2 μs            | 1.2 μs            |
| n = 64                    | 83 ns         | 1.2 μs            | 1.2 μs            |
| n = 128                   | 93 ns         | 1.4 μs            | 1.4 μs            |
| n = 256                   | 111 ns        | 2.0 μs            | 2.0 μs            |
| n = 512                   | 148 ns        | 3.3 μs            | 3.3 μs            |
| n = 1024                  | 228 ns        | 5.9 μs            | 5.9 μs            |
| n = 2048                  | 426 ns        | 11.4 μs           | 11.4 μs           |

| GPU (parallel over SIMDs) | Latency (1x1) | Latency (8x8, LL) | Latency (8x8, RL) |
| ------------------------- | ------------- | ----------------- | ----------------- |
| n = 32                    | 127 ns        | 7.9 μs            | 3.1 μs            |
| n = 64                    | 127 ns        | 7.9 μs            | 3.1 μs            |
| n = 128                   | 127 ns        | 7.9 μs            | 3.1 μs            |
| n = 256                   | 127 ns        | 7.9 μs            | 3.1 μs            |
| n = 512                   | 127 ns        | 7.9 μs            | 3.1 μs            |
| n = 1024                  | 136 ns        | 8.4 μs            | 3.6 μs            |
| n = 2048                  | 154 ns        | 9.3 μs            | 4.5 μs            |

The results are best summarized in the next two tables. The GPU seems to have an advantage for extremely large problem sizes. The latency asymptotically scales with $O((\log{n})^2)$, not $O(n)$.

| Best Latency (8x8) | CPU       | CPU (AMX) | GPU (serial) | GPU (parallel) |
| ------------------ | --------- | --------- | ------------ | -------------- |
| n = 32             | 550 ns    | 530 ns    | 1.2 μs       | 3.1 μs         |
| n = 64             | 620 ns    | 620 ns    | 1.2 μs       | 3.1 μs         |
| n = 128            | 840 ns    | 720 ns    | 1.4 μs       | 3.1 μs         |
| n = 256            | 1.5 μs    | 1.3 μs    | 2.0 μs       | 3.1 μs         |
| n = 512            | 2.7 μs    | 2.3 μs    | 3.3 μs       | 3.1 μs         |
| n = 1024           | 5.1 μs    | 4.3 μs    | 5.9 μs       | 3.6 μs         |
| n = 2048           | 9.7 μs    | 8.2 μs    | 11.4 μs      | 4.5 μs         |

| n   | Best Latency | Processor |
| --- | ------------ | --------- |
| 32  | 530 ns       | CPU       |
| 64  | 620 ns       | CPU       |
| 128 | 720 ns       | CPU       |
| 256 | 1.3 μs       | CPU       |
| 512 | 2.3 μs       | CPU       |
| 1024 | 3.6 μs      | GPU       |
| 2048 | 4.5 μs      | GPU       |

### 1-2-8, 1-2-4-8 Blocking Scheme

Determine how much speedup can be achieved on single-core CPU (Firestorm), with a logarithmically scaling blocking hierarchy. The results are especially important for SBR, which employs a large number of tiny QR factorizations.

This may require an IRL performance benchmark.

## Current Plan

It seems like the bottleneck in bulge chasing could be fixed with AMX GEMV. This should be investigated before trying multi-stage successive band reduction. It is also questionable whether accelerating the panel factorizations is worthwhile at the moment. It could be an important optimization for medium-small semiempirical simulations. However, this part and the bulge chasing are $O(n^2)$.
