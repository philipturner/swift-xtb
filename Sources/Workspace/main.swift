//
//  main.swift
//
//
//  Created by Philip Turner on 5/25/24.
//

import Meshing

// Create an ansatz for the nitrogen atom.
var ansatzDesc = AnsatzDescriptor()
ansatzDesc.atomicNumber = 7
ansatzDesc.netSpin = 1.5
ansatzDesc.position = .zero
ansatzDesc.sizeExponent = 5
let ansatz = Ansatz(descriptor: ansatzDesc)
print(ansatz.orbitals.reduce(0) { $0 + $1.octree.nodes.count })

// Create a mesh from the ansatz.
var meshDesc = MeshDescriptor()
meshDesc.octrees = ansatz.orbitals.map(\.octree)
meshDesc.positions = Array(repeating: .zero, count: ansatz.orbitals.count)
meshDesc.sizeExponent = 3
let mesh = Mesh(descriptor: meshDesc)
