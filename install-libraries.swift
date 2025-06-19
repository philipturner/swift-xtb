import Foundation

// MARK: - Parse Arguments

// Fetch the command-line arguments.
var arguments = CommandLine.arguments

// Parse the script name.
guard arguments.count >= 1 else {
  fatalError("Not enough arguments.")
}
let scriptName = arguments[0]
guard scriptName == "install-libraries.swift" else {
  fatalError("First argument is not a script file.")
}
arguments.removeFirst()

// Parse the otool output.
guard arguments.count >= 1 else {
  fatalError("Not enough arguments.")
}
let otoolOutput = arguments[0]
guard otoolOutput.count > 30 else {
  fatalError("otool output was too short.")
}
arguments.removeFirst()

// Parse the remaining arguments.
var checkOpenBLAS: Bool = false
var checkAccelerate: Bool = false
var reportOpenBLAS: Bool = false

for argument in arguments {
  switch argument {
  case "--check-openblas":
    checkOpenBLAS = true
  case "--check-accelerate":
    checkAccelerate = true
  case "--report-openblas":
    reportOpenBLAS = true
  default:
    fatalError("Unexpected argument: \(argument)")
  }
}

if checkOpenBLAS && checkAccelerate {
  fatalError("--check-openblas and --check-accelerate are mutually exclusive.")
}

// MARK: - Execute Operations

// Split the otool output into lines.
let otoolLines = otoolOutput.split(separator: "\n").map(String.init)

// Check the number of linear algebra library dependencies.
var openblasCount: Int = .zero
var accelerateCount: Int = .zero
for line in otoolLines {
  let openblasRange = line.range(
    of: "openblas", options: .caseInsensitive)
  let accelerateRange = line.range(
    of: "accelerate", options: .caseInsensitive)
  
  if openblasRange != nil {
    openblasCount += 1
  }
  if accelerateRange != nil {
    accelerateCount += 1
  }
}
if checkOpenBLAS {
  guard openblasCount == 1,
        accelerateCount == 0 else {
    fatalError("Unexpected number of linear algebra library dependencies.")
  }
} else if checkAccelerate {
  guard openblasCount == 0,
        accelerateCount == 1 else {
    fatalError("Unexpected number of linear algebra library dependencies.")
  }
}

// Print the OpenBLAS library address to stdout.
if reportOpenBLAS {
  var openblasLine: String?
  for line in otoolLines {
    let openblasRange = line.range(
      of: "openblas", options: .caseInsensitive)
    if openblasRange != nil {
      openblasLine = line
    }
  }
  guard let openblasLine else {
    fatalError("Could not find OpenBLAS line.")
  }
  
  // Get the first instance of each substring in the line.
  let optRange = openblasLine.range(
    of: "/opt", options: .caseInsensitive)
  let dylibRange = openblasLine.range(
    of: ".dylib", options: .caseInsensitive)
  guard let optRange,
        let dylibRange else {
    fatalError("Could not locate library address.")
  }
  
  guard optRange.upperBound < dylibRange.lowerBound else {
    fatalError("Malformatted library address.")
  }
  
  let startIndex = optRange.lowerBound
  let endIndex = dylibRange.upperBound
  let substringRange = startIndex..<endIndex
  let substring = openblasLine[substringRange]
  print(substring)
}
