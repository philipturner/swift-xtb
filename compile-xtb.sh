mkdir .build
cd .build
rm -rf xtb
git clone https://github.com/grimme-lab/xtb

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
rm -rf libxtb_accelerate.dylib
if [ -f "libxtb_accelerate.dylib" ]; then
  echo "Could not remove existing dylib."
  exit -1
fi

# Copy the library to the package directory.
XTB_DIR="$(pwd)/.build/xtb/build"
cp "$XTB_DIR/libxtb.6.7.1.dylib" libxtb_accelerate.dylib
if [ ! -f "libxtb_accelerate.dylib" ]; then
  echo "Could not copy the fresh dylib."
  exit -1
fi

# Homebrew download:
#
#  libxtb_accelerate.dylib:
#    libxtb_accelerate.dylib (compatibility version 6.0.0, current version 6.0.0)
#    /System/Library/Frameworks/Accelerate.framework/Versions/A/Accelerate (compatibility version 0.0.0, current version 0.0.0)
#    /Users/philipturner/Documents/MolecularRenderer/swift-xtb/.build/mctc-lib/0.3.2_1/lib/libmctc-lib.0.dylib (compatibility version 0.0.0, current version 0.0.0)
#    /opt/homebrew/opt/gcc/lib/gcc/current/libgfortran.5.dylib (compatibility version 6.0.0, current version 6.0.0)
#    /opt/homebrew/opt/gcc/lib/gcc/current/libgomp.1.dylib (compatibility version 2.0.0, current version 2.0.0)
#    /usr/lib/libSystem.B.dylib (compatibility version 1.0.0, current version 1351.0.0)
#
# binary size:
# - libxtb_accelerate.dylib: 5.5 MB
# - libmctc-lib.0.dylib: 304 KB

# Compiling from source:
#
#  libxtb_accelerate.dylib:
#    @rpath/libxtb.6.dylib (compatibility version 6.0.0, current version 6.7.1)
#    /opt/homebrew/opt/gcc/lib/gcc/current/libgomp.1.dylib (compatibility version 2.0.0, current version 2.0.0)
#    /opt/homebrew/opt/gcc/lib/gcc/current/libgfortran.5.dylib (compatibility version 6.0.0, current version 6.0.0)
#    /System/Library/Frameworks/Accelerate.framework/Versions/A/Accelerate (compatibility version 1.0.0, current version 4.0.0)
#    /opt/homebrew/opt/gcc/lib/gcc/current/libquadmath.0.dylib (compatibility version 1.0.0, current version 1.0.0)
#    /usr/lib/libSystem.B.dylib (compatibility version 1.0.0, current version 1351.0.0)
#
# binary size:
# - libxtb_accelerate.dylib: 6.4 MB

install_name_tool -id \
  "libxtb_accelerate.dylib" \
  libxtb_accelerate.dylib

echo ""
echo "This code-sign should report success:"
codesign --verify --verbose libxtb_accelerate.dylib
echo ""
