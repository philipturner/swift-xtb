//
//  xTB_Calculator+Properties.swift
//
//
//  Created by Philip Turner on 5/30/24.
//

extension xTB_Calculator {
  /// Numerical accuracy of calculator (in atomic units).
  ///
  /// The default value is 1. The value may range from 1e3 to 1e-4.
  public var accuracy: Float {
    get {
      state.accuracy
    }
    set {
      state.accuracy = newValue
      setAccuracy(newValue)
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
      setMaximumIterations(newValue)
    }
  }
  
  /// Electronic temperature for level filling (in Kelvin).
  public var electronicTemperature: Float {
    get {
      state.electronicTemperature
    }
    set {
      state.electronicTemperature = newValue
      setElectronicTemperature(newValue)
    }
  }
  
  func setAccuracy(_ accuracy: Float) {
    xtb_setAccuracy(
      environment.pointer,
      self.pointer,
      Double(state.accuracy))
  }
  
  func setMaximumIterations(_ maximumIterations: Int) {
    xtb_setMaxIter(
      environment.pointer,
      self.pointer,
      Int32(maximumIterations))
  }
  
  func setElectronicTemperature(_ electronicTemperature: Float) {
    xtb_setElectronicTemp(
      environment.pointer,
      self.pointer,
      Double(electronicTemperature))
  }
}
