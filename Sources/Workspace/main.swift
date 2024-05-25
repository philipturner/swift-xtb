//
//  main.swift
//
//
//  Created by Philip Turner on 5/25/24.
//

import Mechanosynthesis
import Numerics

var ansatzDesc = AnsatzDescriptor()
ansatzDesc.atomicNumber = 1
ansatzDesc.netSpin = 0.5
ansatzDesc.position = .zero
ansatzDesc.sizeExponent = 4
let ansatz = Ansatz(descriptor: ansatzDesc)
print(ansatz.orbitals.count)
print(ansatz.occupations)

let orbital = ansatz.orbitals[0]
print(orbital.octree.nodes.count)

var fragmentCount: Int = .zero
for node in orbital.octree.nodes {
  
}
