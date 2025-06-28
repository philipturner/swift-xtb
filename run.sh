# Prevent crashes with large systems.
export OMP_STACKSIZE="2G"

# Only use the performance cores on macOS.
if [[ "$OSTYPE" == "darwin"* ]]; then
  export OMP_NUM_THREADS=$(sysctl -n hw.perflevel0.physicalcpu)
fi

# Tell the linker where the xTB library is.
export XTB_LIBRARY_PATH="$(pwd)"

# Run in release mode with incremental compilation.
swift run -Xswiftc -Ounchecked
