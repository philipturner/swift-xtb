mkdir .build
cd .build
rm -rf xtb
git clone --single-branch --branch fix-gfnff-output https://github.com/philipturner/xtb

# Flags to explore:
# impact of '-fopenmp'
# try changing to CMake Debug mode
# restrict to 1 OpenMP thread to isolate changes to single-threaded performance
# impact of 'mtune' (https://stackoverflow.com/a/75739441)
# impact of flags to enable loop unrolling/vectorization

cd xtb
cmake -B build \
 -DCMAKE_BUILD_TYPE=Release \
 -DCMAKE_C_COMPILER="gcc-15" \
 -DCMAKE_Fortran_COMPILER="gfortran-15" \
 -DWITH_TBLITE=OFF \
 -DWITH_CPCMX=OFF \
 -DBLA_VENDOR=Apple
make -C build -j8

cd ../ # balance 'cd xtb'
cd ../ # balance 'cd .build'

# Purge the existing dylib to avoid "Permission denied" errors.
rm -rf libxtb.dylib
if [ -f "libxtb.dylib" ]; then
  echo "Could not remove existing dylib."
  exit -1
fi

# Copy the library to the package directory.
XTB_DIR="$(pwd)/.build/xtb/build"
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
