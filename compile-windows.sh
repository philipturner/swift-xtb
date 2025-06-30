# Enter the build folder.
mkdir .build
cd .build



# === Prepare MSYS2 Compiler Environment ===

# # Isolate the process of installing MSYS2.
# rm -rf msys64
# rm -rf msys2-installer.sfx.exe
# curl -L -o "msys2-installer.sfx.exe" "https://github.com/msys2/msys2-installer/releases/download/2025-06-22/msys2-base-x86_64-20250622.sfx.exe"
# ./msys2-installer.sfx.exe

# # Make the current working directory 'msys64' while running these commands.
# # They are a one-time setup procedure before MSYS2 can work properly.
# cd msys64
# usr/bin/bash -leo pipefail %*
# usr/bin/sed -i "s/^CheckSpace/#CheckSpace/g" "/etc/pacman.conf"
# cd ../ # balance 'cd msys64'

# # Isolate the process of installing GNU/CMake/Ninja/openblas.
# msys64/usr/bin/pacman --noconfirm -S "mingw-w64-x86_64-gcc-fortran"
# msys64/usr/bin/pacman --noconfirm -S "mingw-w64-x86_64-cmake"
# msys64/usr/bin/pacman --noconfirm -S "mingw-w64-x86_64-ninja"
# msys64/usr/bin/pacman --noconfirm -S "mingw-w64-x86_64-openblas"



# === Build grimme-lab/xtb ===

# The compiler is installed here.
MINGW64_DIR="$(pwd)/msys64/mingw64"

# Compile the code from source.
# rm -rf xtb
# git clone --single-branch --branch fix-gfnff-output https://github.com/philipturner/xtb
cd xtb
# PATH="$MINGW64_DIR/bin:$PATH" cmake \
#  -B build \
#  -DCMAKE_BUILD_TYPE=Release \
#  -DCMAKE_C_COMPILER="gcc" \
#  -DCMAKE_Fortran_COMPILER="gfortran" \
#  -DWITH_TBLITE=OFF \
#  -DWITH_CPCMX=OFF
# PATH="$MINGW64_DIR/bin:$PATH" ninja \
#   -C build -j4
# PATH="$MINGW64_DIR/bin:$PATH" ninja \
#   -C build test

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

cd ../ # balance 'cd build'
cd ../ # balance 'cd xtb'
cd ../ # balance 'cd .build'



# === Install into directory visible to Swift ===

# The xTB binaries are located here.
XTB_DIR="$(pwd)/.build/xtb/build"
ls "$XTB_DIR"

# Copy dynamic libraries into the folder.
# cp "$MINGW64_DIR/bin/libgfortran-5.dll" "libgfortran-5.dll"
# cp "$XTB_DIR/libxtb.dll" "xtb.dll"


# Copy static libraries into the folder.
cp "$MINGW64_DIR/lib/libgcc_s.a" "libgcc_s_seh-1.lib"
cp "$MINGW64_DIR/lib/libgfortran.a" "libgfortran-5.lib"
cp "$MINGW64_DIR/lib/libgomp.a" "libgomp-1.lib"
cp "$MINGW64_DIR/lib/libopenblas.a" "libopenblas.lib"
cp "$MINGW64_DIR/lib/libwinpthread.a" "libwinpthread-1.lib"
cp "$XTB_DIR/libxtb.a" "xtb.lib"

# Alternative valid command for libgcc_s_seh-1
# cp "$MINGW64_DIR/lib/gcc/x86_64-w64-mingw32/15.1.0/libgcc.a" "libgcc_s_seh-1.lib"
