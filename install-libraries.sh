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

# Iterate over the text line by line

# This is too much of a pain. Use a Swift script for the string processing.
line_id=0
openblas_count=0
accelerate_count=0
IFS=$'\n' # Set Internal Field Separator to newline to handle spaces in lines
while read -r line; do
  line_id=$((line_id + 1))
  
  # Add your desired operations on each 'line' here
  if echo "$line" | grep -iq "openblas"; then
    openblas_count=$((openblas_count + 1))
    
    lines='
    "FreeBSD" {1aac7062-bd59-47ee-9261-2f6aa8d9ef53}
    "Windows 10" {64942de7-beb9-418c-9f52-5befcb6f577b}
    "High Sierra" {07f73e1a-a0c4-4190-ade1-79a2e432b4d6}
    "Test Machine" ee opt9d0953a7-ca2a-4667-8c5b-1a9f550b2956dylib (compatibility)'
    echo "$lines" | awk '{ sub("opt" ,""); sub("dylib", "");  print $3 }'
    
    echo "$line"
    echo -e "deed\ndeath\ndone\n/opt/homebrew/opt/openblas/lib/libopenblas.0.dylib" |awk '/^\/opt.*dylib$/'
    echo -e "deed\ndeath\ndone\n /opt/homebrew/opt/openblas/lib/libopenblas.0.dylib ()" |awk '/\/opt.*dylib/'
    echo -e "      /opt/homebrew/opt/openblas/lib/libopenblas.0.dylib (compatibility version 0.0.0, current version 0.0.0)" |awk '/\/opt.*dylib/'
    echo "$line" |awk '/^opt.*dylib$/'
  fi
  if echo "$line" | grep -iq "accelerate"; then
    accelerate_count=$((accelerate_count + 1))
  fi
done <<< "$otool_output" # Use a "here string" to feed the variable content to the loop
echo ""
echo "Processed $line_id lines."
echo "Found OpenBLAS $openblas_count times."
echo "Found Accelerate $accelerate_count times."
