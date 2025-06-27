//
//  main.swift
//  swift-xtb
//
//  Created by Philip Turner on 6/19/25.
//

import Foundation
import xTB



// GFN-FF
//
// System | Correct Energy  |
// ------ | --------------- |
// 122    |   -57126.394 zJ |
// 222    |  -100762.569 zJ |
// 233    |  -208433.751 zJ |
//
// GFN2-xTB
//
// System | Correct Energy  |
// ------ | --------------- |
// 122    |  -452321.592 zJ |
// 222    |  -841344.658 zJ |
// 233    | -1796496.861 zJ |



let cString1 = getenv("OMP_STACKSIZE")
let cString2 = getenv("OMP_NUM_THREADS")
if let cString1 {
  print("OMP_STACKSIZE:", String(cString: cString1))
}
if let cString2 {
  print("OMP_NUM_THREADS:", String(cString: cString2))
}

// Suppress unwanted output from GFN-FF.
let url = FileManager.default.temporaryDirectory
let path = url.relativePath
let worked = FileManager.default.changeCurrentDirectoryPath(path)
guard worked else {
  fatalError("Could not redirect gfnff_topo directory.")
}
xTB_Environment.verbosity = .muted
//xTB_Environment.setOutput("/dev/null")

// Select the system.
let system: [SIMD4<Float>] = diamondSystem233

// Create the calculator.
var calculatorDesc = xTB_CalculatorDescriptor()
calculatorDesc.atomicNumbers = system.map { UInt8($0.w) }
calculatorDesc.positions = system.map {
  SIMD3($0.x, $0.y, $0.z)
}
calculatorDesc.hamiltonian = .forceField
let calculator = xTB_Calculator(descriptor: calculatorDesc)

// Run just one loop iteration.
var minLatency: Double = 1_000_000
for _ in 0..<10 {
  calculator.molecule.positions = system.map {
    SIMD3($0.x, $0.y, $0.z)
  }
  
  let checkpoint0 = Date()
  let energy = calculator.energy
  let checkpoint1 = Date()
  let latency = checkpoint1.timeIntervalSince(checkpoint0)
  minLatency = min(latency, minLatency)
  
  let formattedLatency = String(format: "%.1f", latency * 1000)
  print()
  print("actual latency:", formattedLatency, "ms")
  
  let formattedEnergy = String(format: "%.3f", energy)
  print("energy:", formattedEnergy, "zJ")
  print("orbitals:", calculator.orbitals.count, calculator.orbitals.eigenvalues.count)
  print("bond orders:", calculator.molecule.bondOrders.count)
  print("charges:", calculator.molecule.charges.count)
  print("forces:", calculator.molecule.forces.count)
  
//  let charges = calculator.molecule.charges
//  print()
//  for charge in charges {
//    print(charge)
//  }
//  print()
  
  print()
  print(xTB_Environment.status)
  print(xTB_Environment.flushErrorStack())
  print(xTB_Environment.status)
  xTB_Environment.show()
  print(xTB_Environment.status)
  print()
}

// Summarize the results
let formattedLatency = String(format: "%.1f", minLatency * 1000)
print()
print("minimum latency:", formattedLatency, "ms")



// Largest system - OpenBLAS

/*
 SCC (total)                   0 d,  0 h,  0 min,  3.643 sec
 SCC setup                      ...        0 min,  0.001 sec (  0.032%)
 Dispersion                     ...        0 min,  0.006 sec (  0.176%)
 classical contributions        ...        0 min,  0.001 sec (  0.023%)
 integral evaluation            ...        0 min,  0.025 sec (  0.686%)
 iterations                     ...        0 min,  3.339 sec ( 91.666%)
 molecular gradient             ...        0 min,  0.263 sec (  7.219%)
 printout                       ...        0 min,  0.007 sec (  0.196%)

         :::::::::::::::::::::::::::::::::::::::::::::::::::::
         ::                     SUMMARY                     ::
         :::::::::::::::::::::::::::::::::::::::::::::::::::::
         :: total energy            -412.064362227315 Eh    ::
         :: gradient norm              0.725113756878 Eh/a0 ::
         :: HOMO-LUMO gap              0.884850477748 eV    ::
         ::.................................................::
         :: SCC energy              -417.574064126950 Eh    ::
         :: -> isotropic ES            0.079789399513 Eh    ::
         :: -> anisotropic ES          0.087013227315 Eh    ::
         :: -> anisotropic XC          0.112741599146 Eh    ::
         :: -> dispersion             -0.739231070380 Eh    ::
         :: repulsion energy           5.509701899634 Eh    ::
         :: add. restraining           0.000000000000 Eh    ::
         :: total charge               0.000000000098 e     ::
         :::::::::::::::::::::::::::::::::::::::::::::::::::::


actual latency: 3647.7 ms
energy: -1796496.862 zJ
orbitals: 776 776

E+G (total)                   0 d,  0 h,  0 min,  0.009 sec
distance/D3 list               ...        0 min,  0.000 sec (  0.765%)
non bonded repulsion           ...        0 min,  0.001 sec (  8.518%)
dCN                            ...        0 min,  0.001 sec (  7.975%)
EEQ energy and q               ...        0 min,  0.001 sec ( 10.861%)
D3                             ...        0 min,  0.004 sec ( 43.356%)
EEQ gradient                   ...        0 min,  0.000 sec (  1.512%)
bonds                          ...        0 min,  0.001 sec ( 11.427%)
bend and torsion               ...        0 min,  0.000 sec (  2.155%)
bonded ATM                     ...        0 min,  0.000 sec (  2.719%)
HB/XB (incl list setup)        ...        0 min,  0.001 sec ( 10.509%)


actual latency: 11.0 ms
energy: -208177.354 zJ
orbitals: 0 0
 */



// Largest System - Accelerate

/*
 SCC (total)                   0 d,  0 h,  0 min,  1.817 sec
 SCC setup                      ...        0 min,  0.001 sec (  0.056%)
 Dispersion                     ...        0 min,  0.007 sec (  0.376%)
 classical contributions        ...        0 min,  0.001 sec (  0.050%)
 integral evaluation            ...        0 min,  0.025 sec (  1.399%)
 iterations                     ...        0 min,  1.517 sec ( 83.510%)
 molecular gradient             ...        0 min,  0.261 sec ( 14.342%)
 printout                       ...        0 min,  0.005 sec (  0.264%)

         :::::::::::::::::::::::::::::::::::::::::::::::::::::
         ::                     SUMMARY                     ::
         :::::::::::::::::::::::::::::::::::::::::::::::::::::
         :: total energy            -412.064362227315 Eh    ::
         :: gradient norm              0.725113756878 Eh/a0 ::
         :: HOMO-LUMO gap              0.884850477748 eV    ::
         ::.................................................::
         :: SCC energy              -417.574064126949 Eh    ::
         :: -> isotropic ES            0.079789399513 Eh    ::
         :: -> anisotropic ES          0.087013227315 Eh    ::
         :: -> anisotropic XC          0.112741599146 Eh    ::
         :: -> dispersion             -0.739231070380 Eh    ::
         :: repulsion energy           5.509701899634 Eh    ::
         :: add. restraining           0.000000000000 Eh    ::
         :: total charge               0.000000000098 e     ::
         :::::::::::::::::::::::::::::::::::::::::::::::::::::


actual latency: 1821.9 ms
energy: -1796496.862 zJ
orbitals: 776 776

E+G (total)                   0 d,  0 h,  0 min,  0.008 sec
distance/D3 list               ...        0 min,  0.000 sec (  0.876%)
non bonded repulsion           ...        0 min,  0.001 sec (  9.245%)
dCN                            ...        0 min,  0.001 sec (  8.467%)
EEQ energy and q               ...        0 min,  0.001 sec (  8.053%)
D3                             ...        0 min,  0.003 sec ( 40.246%)
EEQ gradient                   ...        0 min,  0.000 sec (  1.981%)
bonds                          ...        0 min,  0.001 sec ( 13.464%)
bend and torsion               ...        0 min,  0.000 sec (  2.665%)
bonded ATM                     ...        0 min,  0.000 sec (  2.993%)
HB/XB (incl list setup)        ...        0 min,  0.001 sec ( 11.875%)


actual latency: 9.1 ms
energy: -208177.354 zJ
orbitals: 0 0
 */



// Largest System - Main Branch

/*
 SCC (total)                   0 d,  0 h,  0 min,  1.485 sec
 SCC setup                      ...        0 min,  0.001 sec (  0.064%)
 Dispersion                     ...        0 min,  0.006 sec (  0.387%)
 classical contributions        ...        0 min,  0.000 sec (  0.025%)
 integral evaluation            ...        0 min,  0.023 sec (  1.561%)
 iterations                     ...        0 min,  1.348 sec ( 90.726%)
 molecular gradient             ...        0 min,  0.104 sec (  7.033%)
 printout                       ...        0 min,  0.003 sec (  0.201%)

         :::::::::::::::::::::::::::::::::::::::::::::::::::::
         ::                     SUMMARY                     ::
         :::::::::::::::::::::::::::::::::::::::::::::::::::::
         :: total energy            -412.064362216177 Eh    ::
         :: gradient norm              0.724468099442 Eh/a0 ::
         :: HOMO-LUMO gap              0.884850544720 eV    ::
         ::.................................................::
         :: SCC energy              -417.574064115811 Eh    ::
         :: -> isotropic ES            0.079789399330 Eh    ::
         :: -> anisotropic ES          0.087013227134 Eh    ::
         :: -> anisotropic XC          0.112741599357 Eh    ::
         :: -> dispersion             -0.739231070457 Eh    ::
         :: repulsion energy           5.509701899634 Eh    ::
         :: add. restraining           0.000000000000 Eh    ::
         :: total charge               0.000000000098 e     ::
         :::::::::::::::::::::::::::::::::::::::::::::::::::::


actual latency: 1490.0 ms
energy: -1796496.861 zJ
orbitals: 776 776
 
E+G (total)                   0 d,  0 h,  0 min,  0.007 sec
distance/D3 list               ...        0 min,  0.000 sec (  1.811%)
non bonded repulsion           ...        0 min,  0.000 sec (  4.569%)
dCN                            ...        0 min,  0.001 sec (  9.661%)
EEQ energy and q               ...        0 min,  0.001 sec ( 10.600%)
D3                             ...        0 min,  0.003 sec ( 44.505%)
EEQ gradient                   ...        0 min,  0.000 sec (  2.007%)
bonds                          ...        0 min,  0.001 sec (  8.373%)
bend and torsion               ...        0 min,  0.000 sec (  2.370%)
bonded ATM                     ...        0 min,  0.000 sec (  3.424%)
HB/XB (incl list setup)        ...        0 min,  0.001 sec ( 12.387%)


actual latency: 7.6 ms
energy: -208177.354 zJ
orbitals: 0 0
 */
