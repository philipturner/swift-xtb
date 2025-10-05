# Prevent crashes with large systems.
export OMP_STACKSIZE="2G"

# Only use the performance cores on macOS.
if [[ "$OSTYPE" == "darwin"* ]]; then
  export OMP_NUM_THREADS=$(sysctl -n hw.perflevel0.physicalcpu)
fi

# Run in release mode with incremental compilation.
swift test -Xswiftc -Ounchecked
