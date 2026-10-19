mkdir .build
cd .build
rm -rf xtb
git clone --single-branch --branch fix-gfnff-output https://github.com/philipturner/xtb

# Download the dependency versions that this xtb fork expects. Otherwise, CMake
# may find incompatible versions installed elsewhere (e.g. by conda). And when
# CMake fetches them itself, it ignores the version and takes the latest commit.
rm -rf xtb-deps
mkdir xtb-deps
git clone --depth 1 --branch v0.3.2 https://github.com/grimme-lab/mctc-lib xtb-deps/mctc-lib
git clone --depth 1 --branch v0.5.0 https://github.com/fortran-lang/test-drive xtb-deps/test-drive
XTB_DEPS_DIR="$(pwd)/xtb-deps"

cd xtb
cmake -B build \
 -DCMAKE_BUILD_TYPE=Release \
 -DCMAKE_C_COMPILER="gcc-16" \
 -DCMAKE_Fortran_COMPILER="gfortran-16" \
 -DWITH_TBLITE=OFF \
 -DWITH_CPCMX=OFF \
 -DBLA_VENDOR=Apple \
 -DMCTCLIB_FIND_METHOD=fetch \
 -DTEST_DRIVE_FIND_METHOD=fetch \
 -DFETCHCONTENT_SOURCE_DIR_MCTC-LIB="$XTB_DEPS_DIR/mctc-lib" \
 -DFETCHCONTENT_SOURCE_DIR_TEST-DRIVE="$XTB_DEPS_DIR/test-drive"

# Only build the library. The C API test in this fork fails to link.
make -C build -j8 lib-xtb-shared

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
