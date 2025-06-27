//
//  xTB_CalculatorStorage.swift
//  swift-xtb
//
//  Created by Philip Turner on 5/30/24.
//

import C_xTB

/// All of the state variables inside xTB that must be monitored during a
/// calculation.
struct xTB_CalculatorStorage {
  // Immediately synchronized properties.
  var accuracy: Float = 1.0
  var electronicTemperature: Float = 300
  var maximumIterations: Int = 250
  
  // Lazily synchronized properties.
  var molecule: xTB_Molecule!
  var orbitals: xTB_Orbitals!
}

extension xTB_Calculator {
  public var molecule: xTB_Molecule {
    _read {
      yield storage.molecule!
    }
    _modify {
      yield &storage.molecule!
      invalidateSinglepoint()
    }
  }
  
  public var orbitals: xTB_Orbitals { storage.orbitals! }
  
  /// Numerical accuracy of calculator.
  ///
  /// The default value is 1. The value may range from 1e3 to 1e-4.
  public var accuracy: Float {
    get {
      storage.accuracy
    }
    set {
      storage.accuracy = newValue
      xtb_setAccuracy(
        xTB_Environment.tEnvironment, tCalculator, Double(newValue))
      invalidateSinglepoint()
    }
  }
  
  /// Maximum number of self-consistency iterations.
  ///
  /// The default value is 250.
  ///
  /// > Note: Not available for GFN-FF.
  public var maximumIterations: Int {
    get {
      storage.maximumIterations
    }
    set {
      storage.maximumIterations = newValue
      xtb_setMaxIter(
        xTB_Environment.tEnvironment, tCalculator, Int32(newValue))
      invalidateSinglepoint()
    }
  }
  
  /// Electronic temperature for level filling (in Kelvin).
  ///
  /// The default value is 300 K.
  ///
  /// > Note: Not available for GFN-FF.
  public var electronicTemperature: Float {
    get {
      storage.electronicTemperature
    }
    set {
      storage.electronicTemperature = newValue
      xtb_setElectronicTemp(
        xTB_Environment.tEnvironment, tCalculator, Double(newValue))
      invalidateSinglepoint()
    }
  }
  
  /// Potential energy (in zeptojoules).
  public var energy: Double {
    ensureEnergyCached()
    return results.energy!
  }
}
