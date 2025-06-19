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

// No flags are accepted right now.
guard arguments.count == 0 else {
  fatalError("Got unexpected arguments.")
}

// Objectives:
// - Split the 'otool' string into lines
// - In each line, search for the presence of OpenBLAS and Accelerate
// - Verify the expected number of 'openblas' or 'accelerate' instances
// - Extract the location of 'openblas'
// - Notify the calling program of failures, perhaps through exit code
