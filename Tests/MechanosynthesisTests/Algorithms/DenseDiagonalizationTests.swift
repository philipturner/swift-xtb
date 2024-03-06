import XCTest
import Accelerate // Gate out once there are Swift kernels for CPUs without AMX.
import Numerics
import QuartzCore

final class DenseDiagonalizationTests: XCTestCase {
  // TODO: Create a custom direct diagonalization kernel. Start with the
  // simplest one known, then move on to more advanced solvers.
}
