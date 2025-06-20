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

# Download each binary dependency from the Internet.
curl -OsL "https://github.com/grimme-lab/homebrew-qc/releases/download/mctc-lib-0.3.2_1/mctc-lib-0.3.2_1.arm64_sequoia.bottle.tar.gz"
curl -OsL "https://github.com/grimme-lab/homebrew-qc/releases/download/xtb-6.7.1/xtb-6.7.1.arm64_sequoia.bottle.tar.gz"
tar -xzf "mctc-lib-0.3.2_1.arm64_sequoia.bottle.tar.gz"
tar -xzf "xtb-6.7.1.arm64_sequoia.bottle.tar.gz"

# ## How the binaries should look:
#
#  libxtb.6.dylib:
#    /opt/homebrew/opt/xtb/lib/libxtb.6.dylib (compatibility version 6.0.0, current version 6.0.0)
#    /opt/homebrew/opt/openblas/lib/libopenblas.0.dylib (compatibility version 0.0.0, current version 0.0.0)
#    /opt/homebrew/opt/mctc-lib/lib/libmctc-lib.0.dylib (compatibility version 0.0.0, current version 0.0.0)
#    /opt/homebrew/opt/gcc/lib/gcc/current/libgfortran.5.dylib (compatibility version 6.0.0, current version 6.0.0)
#    /opt/homebrew/opt/gcc/lib/gcc/current/libgomp.1.dylib (compatibility version 2.0.0, current version 2.0.0)
#    /usr/lib/libSystem.B.dylib (compatibility version 1.0.0, current version 1351.0.0)
#
#  libmctc-lib.0.dylib:
#    /opt/homebrew/opt/mctc-lib/lib/libmctc-lib.0.dylib (compatibility version 0.0.0, current version 0.0.0)
#    /opt/homebrew/opt/gcc/lib/gcc/current/libgfortran.5.dylib (compatibility version 6.0.0, current version 6.0.0)
#    /usr/lib/libSystem.B.dylib (compatibility version 1.0.0, current version 1351.0.0)
#
#  xtb:
#    /opt/homebrew/opt/openblas/lib/libopenblas.0.dylib (compatibility version 0.0.0, current version 0.0.0)
#    /opt/homebrew/opt/mctc-lib/lib/libmctc-lib.0.dylib (compatibility version 0.0.0, current version 0.0.0)
#    /opt/homebrew/opt/gcc/lib/gcc/current/libgfortran.5.dylib (compatibility version 6.0.0, current version 6.0.0)
#    /opt/homebrew/opt/gcc/lib/gcc/current/libgomp.1.dylib (compatibility version 2.0.0, current version 2.0.0)
#    /usr/lib/libSystem.B.dylib (compatibility version 1.0.0, current version 1351.0.0)

LIBXTB_PATH="xtb/6.7.1/lib/libxtb.6.dylib"
libxtb_output=$(otool -L "$LIBXTB_PATH")
echo "$libxtb_output"
#swift "../compile-libraries.swift" "$libxtb_output" openblas

LIBMCTC_LIB_PATH="mctc-lib/0.3.2_1/lib/libmctc-lib.0.dylib"
libmctc_lib_output=$(otool -L "$LIBMCTC_LIB_PATH")
echo "$libmctc_lib_output"

XTB_PATH="xtb/6.7.1/bin/xtb"
xtb_output=$(otool -L "$XTB_PATH")
echo "$xtb_output"
