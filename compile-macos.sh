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

# The GCC runtime libraries, from the Homebrew installation that compiled xTB.
DEPENDENCIES=(
  "libgcc_s.1.1.dylib"
  "libgfortran.5.dylib"
  "libgomp.1.dylib"
  "libquadmath.0.dylib"
)

# Make a binary load its dependencies from its own folder, and nowhere else.
#
# Every reference to /opt/homebrew becomes '@rpath', and the only search path
# is '@loader_path'. So the binary fails to load unless the dependency dylibs
# are in the same folder.
make_relocatable() {
  local binary="$1"
  chmod u+w "$binary"
  
  # The first line of 'otool -L' is the binary's own install name.
  if [[ "$binary" == *.dylib ]]; then
    install_name_tool -id "@rpath/$(basename "$binary")" "$binary"
  fi
  for dependency in $(otool -L "$binary" | tail -n +2 | awk '{print $1}'); do
    if [[ "$dependency" == /opt/homebrew/* ]]; then
      install_name_tool -change \
        "$dependency" "@rpath/$(basename "$dependency")" "$binary"
    fi
  done
  
  local rpaths=$(otool -l "$binary" \
    | grep -A2 LC_RPATH | grep " path " | awk '{print $2}')
  for rpath in $rpaths; do
    if [ "$rpath" != "@loader_path" ]; then
      install_name_tool -delete_rpath "$rpath" "$binary"
    fi
  done
  if ! echo "$rpaths" | grep -qx "@loader_path"; then
    install_name_tool -add_rpath "@loader_path" "$binary"
  fi
  
  # Editing the load commands invalidates the code signature.
  codesign --force --sign - "$binary" 2> /dev/null
  
  if otool -L "$binary" | grep -q "/opt/homebrew"; then
    echo "$binary still references /opt/homebrew."
    exit -1
  fi
  if otool -l "$binary" | grep -A2 LC_RPATH | grep -q "/opt/homebrew"; then
    echo "$binary still searches /opt/homebrew."
    exit -1
  fi
}

# Package the build folder, like the Windows script. The dependencies sit next
# to the 'xtb' executable.
cd "$XTB_DIR"
for dependency in "${DEPENDENCIES[@]}"; do
  rm -rf "$dependency"
  cp "$(gfortran-16 -print-file-name="$dependency")" "$dependency"
  make_relocatable "$dependency"
done
make_relocatable "xtb"
make_relocatable "libxtb.6.7.1.dylib"

echo ""
echo "Should see a message containing 'University of Bonn':"
./xtb --version
echo ""

cd ../ # balance 'cd build'
cd ../ # balance 'cd xtb'
cd ../ # balance 'cd .build'

# Purge the existing dylibs to avoid "Permission denied" errors.
rm -rf libxtb.dylib
for dependency in "${DEPENDENCIES[@]}"; do
  rm -rf "$dependency"
done
if [ -f "libxtb.dylib" ]; then
  echo "Could not remove existing dylib."
  exit -1
fi

# Copy the libraries to the package directory.
cp "$XTB_DIR/libxtb.6.7.1.dylib" libxtb.dylib
for dependency in "${DEPENDENCIES[@]}"; do
  cp "$XTB_DIR/$dependency" "$dependency"
done
if [ ! -f "libxtb.dylib" ]; then
  echo "Could not copy the fresh dylib."
  exit -1
fi

# The install name must match the file name, for '@rpath' lookups.
install_name_tool -id "@rpath/libxtb.dylib" libxtb.dylib
codesign --force --sign - libxtb.dylib 2> /dev/null

echo ""
echo "This code-sign should report success:"
for binary in libxtb.dylib "${DEPENDENCIES[@]}"; do
  codesign --verify --verbose "$binary"
done
echo ""
