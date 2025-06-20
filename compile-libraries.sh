# Goal:
# - To start off, replicate the functionality of Homebrew.
#   - Download the raw binary files that Grimme Lab is hosting on GitHub.
#   - Fix the linker issues.
#   - Do everything in the '.build' folder of this project. No need to copy to
#     the top-level project folder, except for libraries that the Swift
#     compiler must be able to see / have transparent access to modifications.
#   - Get the 'xtb' executable running from this mutated install.
# - Compile all necessary dependencies from source, using v0.3.2 (before the
#   bug fix I'm seeking).
