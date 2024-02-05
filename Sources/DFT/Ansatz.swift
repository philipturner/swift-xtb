//
//  Ansatz.swift
//
//
//  Created by Philip Turner on 2/4/24.
//

import OpenCL

// For reference: an implementation of the ansatz that reproduces atomic
// radii very well.
// - Don't store a stack for the Laguerre polynomials. Just the last two
//   terms.
// - We'll need to slightly modify it for core electrons. It only handles
//   the valence electrons for now, to reproduce experimental data.
// - Make the bottom term scale with Z_eff = Z_inner + sqrt(Z_valence).
//   Use sqrt-scaling to emulate subshell screening as well. See whether
//   that reproduces ion core radii.
//   - This task can be done without creating the GPU code.
/*
 func createGeometry() -> [Entity] {
   for Z in 1...118 {
     var n: Int
     var Zeff: Int
     if Z <= 2 {
       n = 1
       Zeff = Z
     } else if Z <= 10 {
       n = 2
       Zeff = Z - 2
     } else if Z <= 18 {
       n = 3
       Zeff = Z - 10
     } else if Z <= 36 {
       n = 4
       if Z <= 20 {
         Zeff = Z - 18
       } else if Z <= 30 {
         Zeff = 2
       } else {
         Zeff = 2 + (Z - 30)
       }
     } else if Z <= 54 {
       n = 5
       if Z <= 38 {
         Zeff = Z - 36
       } else if Z <= 48 {
         Zeff = 2
       } else {
         Zeff = 2 + (Z - 48)
       }
     } else if Z <= 86 {
       n = 6
       if Z <= 56 {
         Zeff = Z - 54
       } else if Z <= 80 {
         Zeff = 2
       } else {
         Zeff = 2 + (Z - 80)
       }
     } else {
       n = 7
       if Z <= 88 {
         Zeff = Z - 86
       } else if Z <= 112 {
         Zeff = 2
       } else {
         Zeff = 2 + (Z - 112)
       }
     }
     
     var correction = Float(n) * Float(n) / Float(Zeff)
     correction = correction.squareRoot()
     var normalization = correction * correction * correction
     normalization = normalization.squareRoot()
     
     let numbers = QuantumNumbers(n: n, l: 0, m: 0, Z: Zeff)
     var integral: Double = .zero
     var radius: Double = .zero
     while radius < 300 {
       let dr: Double = 0.1
       defer { radius += dr }
       
       let coordinates = SphericalCoordinates(
         cartesian: SIMD3(Float(radius) * correction, 0, 0))
       let waveFunction = Double(hydrogenWaveFunction(
         numbers: numbers, coordinates: coordinates) * normalization)
       
       let observable: Double = waveFunction * radius * waveFunction
       integral += observable * radius * radius * dr
     }
     integral *= 4 * Double.pi
     print("\(integral),")
   }
   exit(0)
 }

 // MARK: - Hydrogen Wave Function

 func laguerrePolynomial(
   alpha: Float, n: Int
 ) -> (_ x: Float) -> Float {
   if n == 0 {
     return { _ in 1 }
   } else if n == 1 {
     return { x in 1 + alpha - x }
   } else if n >= 1 {
     return { x in
       var stack: [Float] = []
       stack.append(1)
       stack.append(1 + alpha - x)
       
       for k in 1...(n - 1) {
         let coeffLeft = Float(2 * k + 1) + alpha - x
         let coeffRight = -(Float(k) + alpha)
         let numerator = coeffLeft * stack[k] + coeffRight * stack[k - 1]
         let denominator = Float(k + 1)
         stack.append(numerator / denominator)
       }
       return stack.last!
     }
   }
   
   fatalError("Unsupported value for n.")
 }

 func sphericalHarmonic(
   l: Int, m: Int
 ) -> (_ theta: Float, _ phi: Float) -> Float {
   if l == 0 {
     if m == 0 {
       return { _, _ in
         1.0 / 2 * (1 / Float.pi).squareRoot()
       }
     }
   } else if l == 1 {
     if m == 0 {
       return { theta, _ in
         1.0 / 2 * (3 / Float.pi).squareRoot() * Float.cos(theta)
       }
     }
   }
   
   fatalError("Unsupported value for m or l.")
 }

 func factorial(_ x: Int) -> Int {
   guard x >= 0 else {
     fatalError("Cannot take factorial of negative number.")
   }
   if x == 0 {
     return 1
   } else {
     var output = x
     var counter = x - 1
     while counter > 0 {
       output *= counter
       counter -= 1
     }
     return output
   }
 }

 struct QuantumNumbers {
   var n: Int
   var l: Int
   var m: Int
   var Z: Int
 }

 struct SphericalCoordinates {
   var r: Float
   var phi: Float
   var theta: Float
   
   init(cartesian: SIMD3<Float>) {
     r = (cartesian * cartesian).sum().squareRoot()
     if r.magnitude < .leastNormalMagnitude {
       phi = 0
       theta = 0
     } else {
       // in the physics convention, phi is theta and theta is phi
       phi = Float.atan2(y: cartesian.y, x: cartesian.x)
       theta = Float.acos(cartesian.z / r)
     }
   }
 }

 func hydrogenWaveFunction(
   numbers: QuantumNumbers,
   coordinates: SphericalCoordinates
 ) -> Float {
   let R = { (r: Float) -> Float in
     let numerator = factorial(numbers.n - numbers.l - 1)
     let denominator = 2 * numbers.n * factorial(numbers.n + numbers.l)
     var normalizationFactor = Float(numerator) / Float(denominator)
     
     let shellPart = Float(2 * numbers.Z) / Float(numbers.n)
     normalizationFactor *= shellPart * shellPart * shellPart
     normalizationFactor.formSquareRoot()
     
     let shellRadiusPart = shellPart * r
     let L = laguerrePolynomial(
       alpha: Float(2 * numbers.l + 1),
       n: numbers.n - numbers.l - 1)
     
     return normalizationFactor
     * exp(-shellRadiusPart / 2)
     * pow(shellRadiusPart, Float(numbers.l))
     * L(shellRadiusPart)
   }
   
   let Y = sphericalHarmonic(l: numbers.l, m: numbers.m)
   let magnitude = R(coordinates.r) * Y(coordinates.theta, coordinates.phi)
   let parity = pow(-1, Float(numbers.l))
   return parity * magnitude
 }

 */
