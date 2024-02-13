import XCTest
import DFT

final class DFTTests: XCTestCase {
  func testExample() throws {
    runOrthogonalizationExperiment()
  }
}

func runOrthogonalizationExperiment() {
  // Hamiltonian: [-0.5 ∇^2 + (-1) / |x|] Ψ = E Ψ
  // Grid bounds: 0 Bohr - 3 Bohr
  // Grid spacing: 0.03 Bohr
  // 100 cells, 10 electrons
  //
  // Boundary conditions:
  // - wavefunction left of 0 Bohr equals wavefunction immediately right
  // - wavefunction right of 3 Bohr is zero
  
  var waveFunctions: [[Float]] = []
  let gridSpacing: Float = 0.03
  
  srand48(2021)
  for _ in 0..<10 {
    var waveFunction: [Float] = []
    for _ in 0..<100 {
      let randomNumber = Float(drand48())
      waveFunction.append(2 * randomNumber - 1)
    }
    
    var sum: Double = .zero
    for fragment in waveFunction {
      sum += Double(fragment * fragment * gridSpacing)
    }
    let normalizationFactor = Float((1 / sum).squareRoot())
    for cellID in 0..<100 {
      waveFunction[cellID] *= normalizationFactor
    }
    waveFunctions.append(waveFunction)
  }
  
  func formatString(_ value: Double) -> String {
    let repr = String(format: "%.3f", value)
    if !repr.starts(with: "-") {
      return " " + repr
    } else {
      return repr
    }
  }
  
  func reportState() {
    print()
    print(
      "  1   ", "|",
      "  x   ", "|",
      "-∇^2/2", "|",
      "-1/|x|", "|",
      "Ψ H Ψ ")
    for electronID in 0..<10 {
      let waveFunction = waveFunctions[electronID]
      var observable1: Double = .zero
      var observableX: Double = .zero
      var observableKinetic: Double = .zero
      var observablePotential: Double = .zero
      var observableHamiltonian: Double = .zero
      
      for cellID in 0..<100 {
        let value = waveFunction[cellID]
        let x = gridSpacing * (Float(cellID) + 0.5)
        let leftValue = (cellID > 0) ? waveFunction[cellID - 1] : value
        let rightValue = (cellID < 99) ? waveFunction[cellID + 1] : 0
        
        let derivativeLeft = (value - leftValue) / gridSpacing
        let derivativeRight = (rightValue - value) / gridSpacing
        let laplacian = (derivativeRight - derivativeLeft) / gridSpacing
        let corePotential = -1 / x.magnitude
        
        let d3r = gridSpacing
        observable1 += Double(value * 1 * value * d3r)
        observableX += Double(value * x * value * d3r)
        observableKinetic += Double(value * -0.5 * laplacian * d3r)
        observablePotential += Double(value * corePotential * value * d3r)
        
        let hamiltonian = -0.5 * laplacian + corePotential * value
        observableHamiltonian += Double(value * hamiltonian * d3r)
      }
      print(
        formatString(observable1), "|",
        formatString(observableX), "|",
        formatString(observableKinetic), "|",
        formatString(observablePotential), "|",
        formatString(observableHamiltonian))
    }
    
    print()
    for cellID in 0..<100 {
      var output: String = ""
      for electronID in 0..<10 {
        let value = waveFunctions[electronID][cellID]
        output += formatString(Double(value)) + " "
      }
      print(output)
    }
  }
  
  reportState()
  
  var newWaveFunctions: [[Float]] = []
  for electronID in 0..<10 {
    let Ψ = waveFunctions[electronID]
    var HΨ: [Float] = []
  }
  waveFunctions = newWaveFunctions
}
