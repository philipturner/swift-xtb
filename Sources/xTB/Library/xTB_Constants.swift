//
//  xTB_Constants.swift
//
//
//  Created by Philip Turner on 5/29/24.
//

// This file contains constants for converting between nanomechanical (SI)
// units and atomic units. The public API converts all quantities to SI for
// ease of use. However, it should be possible to back-convert to atomic units
// when comparing against DFT literature data. That is the reason these
// constants are exposed to the public API.

public let xTB_NmPerBohr: Double = 0.0529177
public let xTB_BohrPerNm: Double = 1 / xTB_NmPerBohr

public let xTB_ZJPerHartree: Double = 4359.7482
public let xTB_HartreePerZJ: Double = 1 / xTB_ZJPerHartree
