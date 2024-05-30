//
//  xTB.swift
//
//
//  Created by Philip Turner on 5/29/24.
//

// Goal: Dedicated Swift wrapper over the C API to xTB, similar to
// swift-openmm. However, the API is designed specifically to facilitate
// nanomechanical engineering. It changes the binary representation and units
// for the exposed properties accordingly.
//
// In addition, reverse-engineer the FIRE minimizer used by the 'xtb' program.
// Reimplement it in Swift and keep the source code handy for future reference.
//
// Finally, the module should define some constants to ease the process of
// converting from nanomechanical to SI units. OpenMM defines some constants
// in the OpenMM C wrapper. Therefore, there is a precedent for API naming
// conventions.
