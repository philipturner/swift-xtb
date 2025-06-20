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
mkdir .build
cd .build

# Links referenced in the Homebrew installer:
#
# https://github.com/grimme-lab/homebrew-qc/releases/download/mctc-lib-0.3.2_1/mctc-lib-0.3.2_1.arm64_sequoia.bottle.tar.gz
# https://github.com/grimme-lab/homebrew-qc/releases/download/multicharge-0.3.0/multicharge-0.3.0.arm64_sequoia.bottle.tar.gz
# https://github.com/grimme-lab/homebrew-qc/releases/download/dftd4-3.7.0/dftd4-3.7.0.arm64_sequoia.bottle.tar.gz
# https://github.com/grimme-lab/homebrew-qc/releases/download/xtb-6.7.1/xtb-6.7.1.arm64_sequoia.bottle.tar.gz

curl -OsL "https://github.com/grimme-lab/homebrew-qc/releases/download/mctc-lib-0.3.2_1/mctc-lib-0.3.2_1.arm64_sequoia.bottle.tar.gz"
curl -OsL "https://github.com/grimme-lab/homebrew-qc/releases/download/xtb-6.7.1/xtb-6.7.1.arm64_sequoia.bottle.tar.gz"
tar -xzf "mctc-lib-0.3.2_1.arm64_sequoia.bottle.tar.gz"
tar -xzf "xtb-6.7.1.arm64_sequoia.bottle.tar.gz"
