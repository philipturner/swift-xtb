//
//  xTB.swift
//
//
//  Created by Philip Turner on 5/29/24.
//

// Goal: Create an API almost identical to MM4ForceField, just
// without requiring an equivalent to 'MM4Parameters'.
//
// Design Specifications:
//
// It should include a built-in energy minimizer:
// - likely implemented on the Swift side
// - copying the FIRE parameters used in the 'xtb' program
//
// It should automatically convert between nanomechanical (SI) units and atomic
// units. The public interface should be in nanomechanical units.
//
// 64-bit quantities should be exposed as 32-bit in the public API.

// First step: Dedicated Swift wrapper over the C API, similar to
// swift-openmm.
//
// The precursor will leave 64-bit quantities unchanged.
