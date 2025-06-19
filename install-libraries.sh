# Step into the build folder, to isolate the immediate effects of file
# processing operations from the top-level folder.
mkdir .build
cd .build

# Download each binary dependency from the Internet.
#
# Installing 'xtb' may take ~30 minutes, if you haven't updated Homebrew
# dependencies in a long time.
brew install xtb
