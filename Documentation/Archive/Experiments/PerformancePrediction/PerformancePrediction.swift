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
