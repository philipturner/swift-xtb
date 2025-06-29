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
#   - libgcc_s_seh-1.dll
#   - KERNEL32.dll
#   - msvcrt.dll
#   - libwinpthread-1.dll
# libopenblas.dll
# - installed at MSYS2 MINGW64
# - dependencies:
#   - libgcc_s_seh-1.dll
#   - libgfortran-5.dll
#   - KERNEL32.dll
#   - msvcrt.dll
#   - libgomp-1.dll
#
# === Secondary dependencies ===
#
# libquadmath-0.dll
# - installed at MSYS2 MINGW64
# - dependencies:
#   - libgcc_s_seh-1.dll
#   - KERNEL32.dll
#   - msvcrt.dll
# ADVAPI32.dll
# - not seen in any folders
# libwinpthread-1.dll
# - installed at both Git Bash and MSYS2 MINGW64
# - dependencies:
#   - KERNEL32.dll
#   - msvcrt.dll
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

# Purge any copied binaries.
rm -rf "libgcc_s_seh-1.dll"
rm -rf "libgfortran-5.dll"
rm -rf "libgomp-1.dll"
rm -rf "libopenblas.dll"
rm -rf "libquadmath-0.dll"
rm -rf "libwinpthread-1.dll"

# Copy the binaries into the folder.
export MINGW64_DIR="/c/msys64/mingw64"
cp "$MINGW64_DIR/bin/libgcc_s_seh-1.dll" "libgcc_s_seh-1.dll"
cp "$MINGW64_DIR/bin/libgfortran-5.dll" "libgfortran-5.dll"
cp "$MINGW64_DIR/bin/libgomp-1.dll" "libgomp-1.dll"
cp "$MINGW64_DIR/bin/libopenblas.dll" "libopenblas.dll"
cp "$MINGW64_DIR/bin/libquadmath-0.dll" "libquadmath-0.dll"
cp "$MINGW64_DIR/bin/libwinpthread-1.dll" "libwinpthread-1.dll"

# PATH for Git Bash terminal:
# /c/Users/phili/bin:/mingw64/bin:/usr/local/bin:/usr/bin:/bin:/mingw64/bin:/usr/bin:/c/Users/phili/bin:/c/Windows/system32:/c/Windows:/c/Windows/System32/Wbem:/c/Windows/System32/WindowsPowerShell/v1.0:/c/Windows/System32/OpenSSH:/cmd:/c/ProgramData/chocolatey/bin:/c/Program Files (x86)/Windows Kits/10/Windows Performance Toolkit:/c/Program Files/PowerShell/7:/c/Program Files (x86)/Windows Kits/8.1/Windows Performance Toolkit:/c/Users/phili/AppData/Local/Microsoft/WindowsApps:/c/Users/phili/miniforge3:/c/Users/phili/miniforge3/python:/c/tools/dart-sdk/bin:/c/Users/phili/AppData/Local/Pub/Cache/bin:/c/Users/phili/AppData/Local/Programs/Microsoft VS Code/bin:/c/Users/phili/AppData/Local/Programs/Swift/Runtimes/6.1.0/usr/bin:/c/Users/phili/AppData/Local/Programs/Swift/Toolchains/6.1.0+Asserts/usr/bin:/usr/bin/vendor_perl:/usr/bin/core_perl

# PATH for MSYS2 MINGW64 terminal:
# /mingw64/bin:/usr/local/bin:/usr/bin:/bin:/c/Windows/System32:/c/Windows:/c/Windows/System32/Wbem:/c/Windows/System32/WindowsPowerShell/v1.0/:/usr/bin/site_perl:/usr/bin/vendor_perl:/usr/bin/core_perl

# One method that works:
# export PATH="/c/msys64/mingw64/bin:$PATH"
# ./xtb.exe --version

# Attempt to isolate the dependent DLLs:
# export MINGW64_DIR="/c/msys64/mingw64"

# rm -rf xtb_exe_dir
# mkdir xtb_exe_dir
# cp xtb.exe xtb_exe_dir/xtb.exe
# cp -r "$MINGW64_DIR/bin" xtb_exe_dir/bin
# export PATH="$(pwd)/xtb_exe_dir/bin:$PATH"

# cd xtb_exe_dir
# ./xtb.exe --version

# Deletable libraries:
# - 2to3 through libgcc_s_seh-1.dll
# - libgmp-10.dll through libgnutlsxx-30.dll
# - libhistory8.dll through libngtcp2-16.dll
# - libp11-kit-0.dll through libpython3.dll
# - libreadline8.dll through libssl-3-x64.dll
# - libstdc++-6.dll through libtermcap-0.dll
# - libtre-5.dll through libuv-1.dll
# - libzstd.dll through zstd

# Required libraries:
# - libgfortran-5.dll
# - libgomp-1.dll
# - libopenblas.dll
# - libquadmath-0.dll
# - libwinpthread-1.dll

./xtb.exe --version
