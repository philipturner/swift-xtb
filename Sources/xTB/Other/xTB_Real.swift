//
//  xTB_Real.swift
//  swift-xtb
//
//  Created by Philip Turner on 10/19/26.
//

protocol xTB_Real: BinaryFloatingPoint, Codable, SIMDScalar { }

extension Float: xTB_Real { }

extension Double: xTB_Real { }
