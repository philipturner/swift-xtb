//
//  xTB_Calculator.swift
//
//
//  Created by Philip Turner on 5/29/24.
//

/// A configuration for a singlepoint calculator.
public struct xTB_CalculatorDescriptor {
  /// Required. The number of protons in each atom's nucleus.
  public var atomicNumbers: [UInt8]?
  
  /// Required. The backing calculation environment.
  public var environment: xTB_Environment?
  
  /// Required. The parametrized method for evaluating forces.
  ///
  /// The default value is GFN2-xTB.
  public var hamiltonian: xTB_Hamiltonian = .tightBinding(2)
  
  /// Required. The net charge of the system.
  ///
  /// The default value is zero.
  public var netCharge: Float = .zero
  
  /// Required. The net spin of the system.
  ///
  /// The default value is zero.
  public var netSpin: Float = .zero
  
  /// Optional. The position of each atom's nucleus (in nm).
  ///
  /// When using GFN-FF, the positions are required to initialize force field
  /// parameters. When using tight binding, positions can be specified
  /// after initialization.
  public var positions: [SIMD3<Float>]?
  
  public init() {
    
  }
}

/// Singlepoint calculator.
public class xTB_Calculator {
  var pointer: xtb_TCalculator
  
  var updateRecord = xTB_UpdateRecord()
  
  var storage: xTB_CalculatorStorage
  
  /// Create new calculator object.
  public init(descriptor: xTB_CalculatorDescriptor) {
    guard let calc = xtb_newCalculator() else {
      fatalError("Could not create new xTB_Calculator.")
    }
    self.pointer = calc
    
    // Create the storage.
    guard let environment = descriptor.environment else {
      fatalError("Descriptor was incomplete.")
    }
    let molecule = xTB_Molecule(descriptor: descriptor)
    storage = xTB_CalculatorStorage(
      environment: environment, molecule: molecule)
    
    // Initialize the parameters.
    if let positions = descriptor.positions {
      self.positions = positions
    }
    flushUpdateRecord()
    loadHamiltonian(descriptor.hamiltonian)
  }
  
  /// Delete calculator object.
  deinit {
    xtb_delCalculator(&pointer)
  }
  
  /// Set numerical accuracy of calculator in the range of 1000 to 0.0001
  public var accuracy: Float {
    get {
      storage.accuracy
    }
    set {
      storage.accuracy = newValue
      updateRecord.accuracy = true
    }
  }
}

extension xTB_Calculator {
  func setAccuracy(_ accuracy: Float) {
    xtb_setAccuracy(
      storage.environment.pointer,
      pointer,
      Double(accuracy))
  }
  
  func setMaximumIterations(_ maximumIterations: Int) {
    xtb_setMaxIter(
      storage.environment.pointer,
      pointer,
      Int32(maximumIterations))
  }
  
  func setElectronicTemperature(_ electronicTemperature: Float) {
    xtb_setElectronicTemp(
      storage.environment.pointer,
      pointer,
      Double(electronicTemperature))
  }
}
