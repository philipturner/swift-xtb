//
//  main.swift
//
//
//  Created by Philip Turner on 5/25/24.
//

import Mechanosynthesis
import Numerics

// Create an ansatz for the nitrogen atom.
var ansatzDesc = AnsatzDescriptor()
ansatzDesc.atomicNumber = 7
ansatzDesc.netSpin = 1.5
ansatzDesc.position = .zero
ansatzDesc.sizeExponent = 5
let ansatz = Ansatz(descriptor: ansatzDesc)

// Create a mesh from the ansatz.
var meshDesc = MeshDescriptor()
meshDesc.octrees = ansatz.orbitals.map(\.octree)
meshDesc.positions = [.zero]
meshDesc.sizeExponent = 3
let mesh = Mesh(descriptor: meshDesc)
