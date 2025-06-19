# Download each binary dependency from the Internet.
if [ ! -d "/opt/homebrew/Cellar/xtb" ]; then
  # Installing 'xtb' may take ~30 minutes, if you haven't updated Homebrew
  # dependencies in a long time.
  brew install xtb
else
  echo "xtb is already installed."
fi

# Bash executes 'while read' in a sub-process, meaning it can't modify
# variables in the calling process. The solution is 'lastpipe'.
#
# Source: https://serverfault.com/a/1039073
shopt -s lastpipe

# Check for multiple installations.
installation_count=0
TARGET_DIR="/opt/homebrew/Cellar/xtb"
find "$TARGET_DIR" -maxdepth 1 -type d -print0 | while IFS= read -r -d $'\0' dir; do
  # Exclude the target directory itself.
  if [ "$dir" != "$TARGET_DIR" ]; then
    installation_count=$((installation_count + 1))
  fi
done
if [ "$installation_count" -eq 1 ]; then
  echo "Identified unique version: $(ls "$TARGET_DIR")"
else
  echo "Multiple versions of xtb are installed. Please uninstall xtb."
  exit -1
fi

# Construct the complete URL of where xtb was downloaded.
XTB_DIR="$TARGET_DIR/$(ls "$TARGET_DIR")"
echo "Installation directory: $XTB_DIR"

# Purge the existing dylib to avoid "Permission denied" errors.
rm -rf libxtb.6.dylib
if [ -f "libxtb.6.dylib" ]; then
  echo "Could not remove existing dylib."
  exit -1
fi

# Copy the library to the package directory.
cp "$XTB_DIR/lib/libxtb.6.dylib" libxtb.6.dylib
if [ ! -f "libxtb.6.dylib" ]; then
  echo "Could not copy the fresh dylib."
  exit -1
fi

# Inspect the dylib's binary dependencies.
otool_output=$(otool -L libxtb.6.dylib)
openblas_address=$(swift "install-libraries.swift" \
  "$otool_output" \
  --check-openblas \
  --report-openblas)
echo "openblas_address = $openblas_address"
echo $openblas_address

# Replace OpenBLAS with Accelerate.
install_name_tool -change \
  "$openblas_address" \
  "/System/Library/Frameworks/Accelerate.framework/Versions/A/Accelerate" \
  libxtb.6.dylib

# Inspect the dylib's binary dependencies.
otool_output=$(otool -L libxtb.6.dylib)
swift "install-libraries.swift" \
  "$otool_output" \
  --check-accelerate

# Source: https://stackoverflow.com/a/2989954
install_name_tool -id \
  "YEET.dylib" \
  libxtb.6.dylib

# Running 'otool' invalidates the code signature. This causes the program to
# crash when loading the dylib through 'dlopen'. The solution is to re-sign
# the dylib with an ad-hoc signature.
#
# Source: https://developer.apple.com/forums/thread/747909
codesign -fs - libxtb.6.dylib

rm -rf YEET.dylib
cp libxtb.6.dylib YEET.dylib
