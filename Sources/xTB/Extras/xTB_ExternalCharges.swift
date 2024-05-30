//
//  xTB_ExternalCharges.swift
//
//
//  Created by Philip Turner on 5/30/24.
//

public struct xTB_ExternalCharge {
  /// The partial charge.
  public var charge: Float = .zero
  
  /// This will be rounded to the nearest integer.
  public var chemicalHardness: Float = .zero
  
  public init() {
    
  }
}

/// A configuration for an external charge potential.
public struct xTB_ExternalPotentialDescriptor {
  public var charges: [xTB_ExternalCharge]?
  
  public var positions: [SIMD3<Float>]?
}

/// External charge potential.
public struct xTB_ExternalPotential {
  public init(descriptor: xTB_ExternalPotentialDescriptor) {
    guard let charges = descriptor.charges,
          let positions = descriptor.positions else {
      fatalError("Descriptor was incomplete.")
    }
    guard charges.count == positions.count else {
      fatalError("Size of charges must match size of positions.")
    }
  }
}
