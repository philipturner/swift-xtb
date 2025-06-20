import Foundation

// MARK: - Parse Arguments

// Fetch the command-line arguments.
var arguments = CommandLine.arguments

// Parse the script name.
guard arguments.count >= 1 else {
  fatalError("Not enough arguments.")
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

// Parse the keyword to search for.
guard arguments.count >= 1 else {
  fatalError("Not enough arguments.")
}
let keyword = arguments[0]
arguments.removeFirst()

// Ensure there are no additional arguments.
guard arguments.count == 0 else {
  fatalError("Too many arguments.")
}

// MARK: - Execute Operations

// Split the otool output into lines.
let otoolLines = otoolOutput.split(separator: "\n").map(String.init)

// Check the number of lines where the keyword appears.
var keywordCount: Int = .zero
var keywordLineID: Int?
for lineID in otoolLines.indices {
  let line = otoolLines[lineID]
  
  let homebrewRange = line.range(
    of: "@@HOMEBREW_PREFIX@@", options: [])
  let keywordRange = line.range(
    of: keyword, options: .caseInsensitive)
  
  if homebrewRange != nil,
     keywordRange != nil {
    keywordCount += 1
    keywordLineID = lineID
  }
}
guard keywordCount == 1,
      let keywordLineID else {
  fatalError("Found the \(keywordCount) instances of \(keyword).")
}

// Print the library address to stdout.
do {
  let keywordLine = otoolLines[keywordLineID]
  
  // Get the first instance of each substring in the line.
  let homebrewRange = keywordLine.range(
    of: "@@HOMEBREW_PREFIX@@", options: [])
  let dylibRange = keywordLine.range(
    of: ".dylib", options: [])
  guard let homebrewRange,
        let dylibRange else {
    fatalError("Could not locate library address.")
  }
  guard homebrewRange.upperBound < dylibRange.lowerBound else {
    fatalError("Malformatted library address.")
  }
  
  // Extract the substring.
  let startIndex = homebrewRange.lowerBound
  let endIndex = dylibRange.upperBound
  let substringRange = startIndex..<endIndex
  let substring = keywordLine[substringRange]
  
  // Report the substring to the calling program.
  print(substring)
}
