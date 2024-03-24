//
//  PerformancePrediction.swift
//
//
//  Created by Philip Turner on 3/21/24.
//

/*
 GFLOPS/k with different block sizes
 
 | Block Size | `sgemm_`   | Theoretical | Agreement  |
 | :--------: | :--------: | :---------: | :--------: |
 | 1x1        | 1.6  -> 12 | 23          | assume 50% |
 | 2x2        | 5.7  -> 23 | 46          | assume 50% |
 | 4x4        | 18.2 -> 23 | 46          | assume 50% |
 | 8x8        | 42.5       | 46          | 92.4%      |
 | 16x16      | 172.4      | 183         | 94.2%      |
 | 32x32      | 566.2      | 733         | 77.2%      |
 | 64x64      | 581.9      | 740         | 78.6%      |
 | 128x128    | 592.0      | 1470        | 40.3%      |
 | 256x256    | 995.4      | 1480        | 67.3%      |
 | 512x512    | 1071.1     | 1480        | 72.4%      |
 | 1024x1024  | 942.3      | 1480        | 63.7%      |
 */

// Predicting the performance of QR factorization.

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
