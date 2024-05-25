//
//  main.swift
//
//
//  Created by Philip Turner on 5/25/24.
//

import Mechanosynthesis
import Numerics

// Tasks:
// - Learn how to set up a mesh for the hydrogen atom.
// - Wrap the grid generation code into a function, store in 'Workspace'.

// Create an ansatz for the electronic structure.
var ansatzDesc = AnsatzDescriptor()
ansatzDesc.atomicNumber = 1
ansatzDesc.netSpin = 0.5
ansatzDesc.position = .zero
ansatzDesc.sizeExponent = 4
let ansatz = Ansatz(descriptor: ansatzDesc)

let orbital = ansatz.orbitals[0]
let mesh = OrbitalMesh(orbital: orbital)

// Fill in the voxels.


