//
//  PerformancePrediction.swift
//
//
//  Created by Philip Turner on 3/21/24.
//

import Foundation

// Predicting the performance of 8x8 tile factorization.

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

// MARK: - Script

let argumentN = CommandLine.arguments[1]
let n: Int = Int(argumentN)!
gpuParallel(n: n, rightLooking: false)
gpuParallel(n: n, rightLooking: true)
