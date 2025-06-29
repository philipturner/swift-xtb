mkdir .build
cd .build
rm -rf xtb
git clone --single-branch --branch fix-gfnff-output https://github.com/philipturner/xtb

# Flags to explore:
# impact of '-fopenmp'
# impact of '-O0, -O1, -O2, -O3'
# impact of 'march/mavx' for certain vector instructions
# impact of flags to enable loop unrolling/vectorization

cd xtb
cmake -B build \
 -DCMAKE_BUILD_TYPE=Release \
 -DCMAKE_C_COMPILER="gcc-15" \
 -DCMAKE_Fortran_COMPILER="gfortran-15" \
 -DCMAKE_Fortran_FLAGS="YEET" \
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
