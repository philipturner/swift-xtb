# Isolate the process of installing MSYS2 and MINGW64
#
# For now, invoke the compiler from a regular command line, using the binaries
# already installed on the PC. After that, we can trace it back to CLI commands
# that download the GNU/CMake/Ninja dependencies.
export MINGW64_DIR="/c/msys64/mingw64"

# Enter the build folder.
mkdir .build
cd .build
rm -rf xtb
git clone --single-branch --branch fix-gfnff-output https://github.com/philipturner/xtb

# Compile the code from source.
cd xtb
PATH="$MINGW64_DIR/bin:$PATH" cmake \
 -B build \
 -DCMAKE_BUILD_TYPE=Release \
 -DCMAKE_C_COMPILER="gcc" \
 -DCMAKE_Fortran_COMPILER="gfortran" \
 -DWITH_TBLITE=OFF \
 -DWITH_CPCMX=OFF
PATH="$MINGW64_DIR/bin:$PATH" ninja \
  -C build -j4
PATH="$MINGW64_DIR/bin:$PATH" ninja \
  -C build test

# Purge any copied binaries.
cd build
rm -rf "libgcc_s_seh-1.dll"
rm -rf "libgfortran-5.dll"
rm -rf "libgomp-1.dll"
rm -rf "libopenblas.dll"
rm -rf "libquadmath-0.dll"
rm -rf "libwinpthread-1.dll"

# Copy the binaries into the folder.
cp "$MINGW64_DIR/bin/libgcc_s_seh-1.dll" "libgcc_s_seh-1.dll"
cp "$MINGW64_DIR/bin/libgfortran-5.dll" "libgfortran-5.dll"
cp "$MINGW64_DIR/bin/libgomp-1.dll" "libgomp-1.dll"
cp "$MINGW64_DIR/bin/libopenblas.dll" "libopenblas.dll"
cp "$MINGW64_DIR/bin/libquadmath-0.dll" "libquadmath-0.dll"
cp "$MINGW64_DIR/bin/libwinpthread-1.dll" "libwinpthread-1.dll"

echo ""
echo "Should see a message containing 'University of Bonn':"
./xtb.exe --version
echo ""
