mkdir .build
cd .build
rm -rf xtb
git clone --single-branch --branch swift-xtb-2026 https://github.com/philipturner/xtb

# Always download mctc-lib at the version xtb expects. Otherwise, CMake may
# find an incompatible version installed elsewhere (e.g. by conda).
cd xtb
cmake -B build \
 -DCMAKE_BUILD_TYPE=Release \
 -DCMAKE_C_COMPILER="gcc-16" \
 -DCMAKE_Fortran_COMPILER="gfortran-16" \
 -DWITH_TBLITE=OFF \
 -DWITH_CPCMX=OFF \
 -DWITH_TESTS=OFF \
 -DBLA_VENDOR=Apple \
 -DMCTCLIB_FIND_METHOD=fetch
make -C build -j8

cd ../ # balance 'cd xtb'
cd ../ # balance 'cd .build'

# Check for the fresh dylib before deleting the existing one.
XTB_DIR="$(pwd)/.build/xtb/build"
if [ ! -f "$XTB_DIR/libxtb.6.7.1.dylib" ]; then
  echo "The xTB build did not produce a dylib."
  exit -1
fi

# Purge the existing dylib to avoid "Permission denied" errors.
rm -rf libxtb.dylib
if [ -f "libxtb.dylib" ]; then
  echo "Could not remove existing dylib."
  exit -1
fi

# Copy the library to the package directory.
cp "$XTB_DIR/libxtb.6.7.1.dylib" libxtb.dylib
if [ ! -f "libxtb.dylib" ]; then
  echo "Could not copy the fresh dylib."
  exit -1
fi

install_name_tool -id \
  "libxtb.dylib" \
  libxtb.dylib

echo ""
echo "This code-sign should report success:"
codesign --verify --verbose libxtb.dylib
echo ""
