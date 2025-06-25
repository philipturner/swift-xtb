mkdir .build
cd .build
rm -rf mctc-lib
rm -rf xtb

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
LIBMCTC_LIB_PATH="mctc-lib/0.3.2_1/lib/libmctc-lib.0.dylib"
XTB_PATH="xtb/6.7.1/bin/xtb"

libxtb_output=$(otool -L "$LIBXTB_PATH")
openblas_address=$(swift "../compile-libraries.swift" "$libxtb_output" openblas)
mctc_lib_address=$(swift "../compile-libraries.swift" "$libxtb_output" mctc-lib)
gfortran_address=$(swift "../compile-libraries.swift" "$libxtb_output" gfortran)
gomp_address=$(swift "../compile-libraries.swift" "$libxtb_output" gomp)

install_name_tool -change \
  "$openblas_address" \
  "/opt/homebrew/opt/openblas/lib/libopenblas.0.dylib" \
  "$LIBXTB_PATH"
install_name_tool -change \
  "$mctc_lib_address" \
  "$(pwd)/$LIBMCTC_LIB_PATH" \
  "$LIBXTB_PATH"
install_name_tool -change \
  "$gfortran_address" \
  "/opt/homebrew/opt/gcc/lib/gcc/current/libgfortran.5.dylib" \
  "$LIBXTB_PATH"
install_name_tool -change \
  "$gomp_address" \
  "/opt/homebrew/opt/gcc/lib/gcc/current/libgomp.1.dylib" \
  "$LIBXTB_PATH"
install_name_tool -id \
  "$(pwd)/$LIBXTB_PATH" \
  "$LIBXTB_PATH"
codesign -fs - "$LIBXTB_PATH"

libmctc_lib_output=$(otool -L "$LIBMCTC_LIB_PATH")
gfortran_address=$(swift "../compile-libraries.swift" "$libmctc_lib_output" gfortran)

install_name_tool -change \
  "$gfortran_address" \
  "/opt/homebrew/opt/gcc/lib/gcc/current/libgfortran.5.dylib" \
  "$LIBMCTC_LIB_PATH"
install_name_tool -id \
  "$(pwd)/$LIBMCTC_LIB_PATH" \
  "$LIBMCTC_LIB_PATH"
codesign -fs - "$LIBMCTC_LIB_PATH"

xtb_output=$(otool -L "$XTB_PATH")
openblas_address=$(swift "../compile-libraries.swift" "$xtb_output" openblas)
mctc_lib_address=$(swift "../compile-libraries.swift" "$xtb_output" mctc-lib)
gfortran_address=$(swift "../compile-libraries.swift" "$xtb_output" gfortran)
gomp_address=$(swift "../compile-libraries.swift" "$xtb_output" gomp)

install_name_tool -change \
  "$openblas_address" \
  "/opt/homebrew/opt/openblas/lib/libopenblas.0.dylib" \
  "$XTB_PATH"
install_name_tool -change \
  "$mctc_lib_address" \
  "$(pwd)/$LIBMCTC_LIB_PATH" \
  "$XTB_PATH"
install_name_tool -change \
  "$gfortran_address" \
  "/opt/homebrew/opt/gcc/lib/gcc/current/libgfortran.5.dylib" \
  "$XTB_PATH"
install_name_tool -change \
  "$gomp_address" \
  "/opt/homebrew/opt/gcc/lib/gcc/current/libgomp.1.dylib" \
  "$XTB_PATH"
codesign -fs - "$XTB_PATH"

echo ""
echo "This code-sign should report success:"
codesign --verify --verbose "$LIBXTB_PATH"
codesign --verify --verbose "$LIBMCTC_LIB_PATH"
codesign --verify --verbose "$XTB_PATH"
echo ""
