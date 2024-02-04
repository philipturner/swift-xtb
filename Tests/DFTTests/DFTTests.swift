import XCTest
import DFT
import OpenCL

final class DFTTests: XCTestCase {
  static let queue: CLCommandQueue = {
    let context = CLContext.default!
    let device = CLDevice.default!
    XCTAssert(device.type == .gpu, "OpenCL device was not a GPU.")
    return CLCommandQueue(context: context, device: device)!
  }()
}
