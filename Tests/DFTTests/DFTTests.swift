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

final class RAMTests: XCTestCase {
  func testGPUAddress() throws {
    // Initialize a RAM with 256 words.
    // +41 ms
    
    let context = CLContext.default!
    let ram = RAM(queue: DFTTests.queue, size: 1024)
    XCTAssertNotEqual(UInt(bitPattern: ram.cpuAddress), 0)
    XCTAssertNotEqual(ram.gpuAddress, 0)
    
    // Issue a GPU kernel to ensure the address is correct.
    // +1 ms
    
    let source = """
    typedef struct {
      float inputA;
      global uint* outputAddress;
    } DataStructure;
    
    // Test something where multiple GPU threads perform slightly different
    // computations.
    kernel void performComputation(global uint* ram,
                                   uint offset) {
      global DataStructure* dataStructures = (global DataStructure*)(ram + offset);
      DataStructure dataStructure = dataStructures[0];
    
      uint threadID = get_global_id(0);
      global float* outputFloat = (global float*)dataStructure.outputAddress;
      outputFloat[threadID] = dataStructure.inputA + float(threadID);
    }
    """
    let program = CLProgram(context: context, source: source)!
    do {
      try program.build(options: "-w")
    } catch {
      let log = program.buildLog(device: context.devices!.first!)
      fatalError("Build error: \(log ?? "n/a")")
    }
    let kernel = program.createKernels()![0]
    
    // Encode the kernel.
    let outputOffset: Int = 128
    let outputPointer = (ram.cpuAddress + outputOffset * 4)
      .assumingMemoryBound(to: Float.self)
    
    struct DataStructure {
      var inputA: Float
      var outputAddress: UInt64
    }
    var ramOffset: Int = 64
    let dataStructurePointer = (ram.cpuAddress + ramOffset * 4)
      .assumingMemoryBound(to: DataStructure.self)
    dataStructurePointer.pointee = DataStructure(
      inputA: 3.5, 
      outputAddress: ram.gpuAddress + 128 * 4)
    
    try! kernel.setArgument(ram.buffer, index: 0)
    try! kernel.setArgument(&ramOffset, index: 1, size: 4)
    let queue = CLCommandQueue(
      context: context, device: context.devices!.first!)!
    try! queue.enqueueKernel(kernel, globalSize: CLNDRange(width: 5))
    try! queue.finish()
    
    for elementID in 0..<5 {
      if elementID < 5 {
        XCTAssertEqual(outputPointer[elementID], 3.5 + Float(elementID))
      } else {
        XCTAssertEqual(outputPointer[elementID], 0)
      }
    }
  }
  
  func testLargeAllocation() throws {
    // + 37 ms/GB
    let ram = RAM(queue: DFTTests.queue, size: 1_000_000_000)
    XCTAssertNotEqual(ram.gpuAddress, 0)
  }
}

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
