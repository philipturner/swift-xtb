import XCTest
import DFT
import OpenCL

final class SelfConsistentFieldTests: XCTestCase {
  func testMultiGridLevel() throws {
    let field = SelfConsistentField()
    field.queue = DFTTests.queue
    field.ram = RAM(queue: field.queue, size: 1024)
    
    var level = MultiGridLevel()
    level.address = 100
    level.spacing = 0.25
    level.origin = SIMD3(-1, -1, -1)
    level.size = SIMD3(1, 1, 1)
     
    let source = """
    typedef struct {
      uint address;
      float spacing;
      int3 origin;
      uint3 size;
    } MultiGridLevel;
    
    kernel void testMultiGridLevel(global uint* RAM,
                                   MultiGridLevel level) {
      float value = level.spacing;
      value += float(level.origin.z);
      value += (float)(level.size.z * 5);
      RAM[level.address] = as_uint(value);
    }
    """
    let program = CLProgram(context: field.queue.context!, source: source)!
    try! program.build()
    let kernel = program.createKernels()![0]
    try! kernel.setArgument(field.ram.buffer, index: 0)
    try! kernel.setArgument(
      &level, index: 1, size: MemoryLayout<MultiGridLevel>.stride)
    try! field.queue.enqueueKernel(kernel, globalSize: CLNDRange(width: 1))
    try! field.queue.finish()
    
    let outputPointer = (field.ram.cpuAddress + 100 * 4)
      .assumingMemoryBound(to: Float.self)
    XCTAssertEqual(outputPointer.pointee, 4.25)
  }
}
