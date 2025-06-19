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
  echo "Multiple version of xtb are installed. Please uninstall xtb."
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
echo "$otool_output"
