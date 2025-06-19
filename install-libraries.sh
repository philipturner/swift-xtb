# Step into the build folder, to isolate the immediate effects of file
# processing operations from the top-level folder.
if [ ! -d ".build" ]; then
  mkdir .build
fi
cd .build

echo "Hello world"
