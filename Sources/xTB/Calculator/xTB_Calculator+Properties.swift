//
//  xTB_Calculator+Properties.swift
//
//
//  Created by Philip Turner on 5/30/24.
//

extension xTB_Calculator {
  public var externalCharges: xTB_ExternalCharges {
    _read {
      yield state.externalCharges!
    }
    _modify {
      yield &state.externalCharges!
      updateRecord.externalCharges = true
      invalidateSinglepoint()
    }
  }
  
  public var molecule: xTB_Molecule {
    _read {
      yield state.molecule!
    }
    _modify {
      yield &state.molecule!
      updateRecord.molecule = true
      invalidateSinglepoint()
    }
  }
  
  public var orbitals: xTB_Orbitals { state.orbitals! }
  
  /// Numerical accuracy of calculator (in atomic units).
  ///
  /// The default value is 1. The value may range from 1e3 to 1e-4.
  public var accuracy: Float {
    get {
      state.accuracy
    }
    set {
      state.accuracy = newValue
      xtb_setAccuracy(
        xTB_Environment._environment, _calculator, Double(newValue))
      invalidateSinglepoint()
    }
  }
  
  /// Maximum number of self-consistency iterations.
  ///
  /// The default value is 250.
  public var maximumIterations: Int {
    get {
      state.maximumIterations
    }
    set {
      state.maximumIterations = newValue
      xtb_setMaxIter(
        xTB_Environment._environment, _calculator, Int32(newValue))
      invalidateSinglepoint()
    }
  }
  
  /// Electronic temperature for level filling (in Kelvin).
  public var electronicTemperature: Float {
    get {
      state.electronicTemperature
    }
    set {
      state.electronicTemperature = newValue
      xtb_setElectronicTemp(
        xTB_Environment._environment, _calculator, Double(newValue))
      invalidateSinglepoint()
    }
  }
}

extension xTB_Calculator {
  /// Potential energy (in zJ).
  public var energy: Double {
    ensureEnergyCached()
    return results.energy!
  }
}
