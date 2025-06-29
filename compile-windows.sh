# === Binary dependencies of xtb.exe/libxtb.dll ===
#
# libgcc_s_seh-1.dll
# - installed at both Git Bash and MSYS2 MINGW64
# - dependencies:
#   - KERNEL32.dll
#   - msvcrt.dll
#   - libwinpthread-1.dll
# libgfortran-5.dll
# - installed at MSYS2 MINGW64
# - dependencies:
#   - libquadmath-0.dll
#   - ADVAPI32.dll
#   - libgcc_s_seh-1.dll
#   - KERNEL32.dll
#   - msvcrt.dll
#   - libwinpthread-1.dll
# KERNEL32.dll
# - not seen in any folders
# msvcrt.dll
# - not seen in any folders
# libgomp-1.dll
# - installed at MSYS2 MINGW64
# - dependencies:
#   - TODO
# libopenblas.dll
# - installed at MSYS2 MINGW64
# - dependencies:
#   - TODO
#
# === Secondary dependencies ===
#
# libquadmath-0.dll
# - installed at MSYS2 MINGW64
# - dependencies:
#   - TODO
# ADVAPI32.dll
# - not seen in any folders
# libwinpthread-1.dll
# - installed at both Git Bash and MSYS2 MINGW64
#
# === Commentary ===
#
# Issue with mingw64:
# - Git Bash terminal: C:\Program Files\Git\mingw64
# - MSYS2 MINGW64 terminal: C:\msys64\mingw64
# - VS Code terminal: not even Bash (cannot do 'export' command)

# Enter the build folder
cd .build
cd xtb
cd build

# Copy the binaries into the folder.
export MINGW64_DIR="/c/msys64/mingw64"
cp "$MINGW64_DIR/bin/libgcc_s_seh-1.dll" "libgcc_s_seh-1.dll"
cp "$MINGW64_DIR/bin/libgfortran-5.dll" "libgfortran-5.dll"
cp "$MINGW64_DIR/bin/libgomp-1.dll" "libgomp-1.dll"
cp "$MINGW64_DIR/bin/libopenblas.dll" "libopenblas.dll"

./xtb.exe --version
