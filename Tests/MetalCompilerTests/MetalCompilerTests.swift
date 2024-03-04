import XCTest
import MetalCompiler

final class MetalCompilerTests: XCTestCase {
  // Test some basic functionality of the Metal compiler.
  func testMTLCompiler() throws {
    let compiler = MTLCompiler()
    compiler.setXcodePath("/Applications/Xcode.app")
    XCTAssertGreaterThan(
      compiler.buildProductsDirectory.relativePath.count, 0)
  }
  
  // Test a basic tutorial that's often used for OpenCL.
  func testVectorAddition() throws {
    let compiler = MTLCompiler()
    let library = compiler.compile("""
      kernel void vectorAddition(
        device float *A [[buffer(0)]],
        device float *B [[buffer(1)]],
        device float *C [[buffer(2)]],
        uint tid [[thread_position_in_grid]]
      ) {
        C[tid] = A[tid] + B[tid];
      }
      """)
    
    let device = MTLCreateSystemDefaultDevice()!
    let function = library.makeFunction(name: "vectorAddition")!
    let pipeline = try! device.makeComputePipelineState(function: function)
    
    let commandQueue = device.makeCommandQueue()!
    let commandBuffer = commandQueue.makeCommandBuffer()!
    let buffer = device.makeBuffer(length: 1024)!
    let contents = buffer.contents().assumingMemoryBound(to: Float.self)
    
    // Define the input operands.
    var operandA: [Float] = [1, 2, 3.3, -9, 1, 5, 59]
    var operandB: [Float] = [-9, 7, 72, 0.123, 5.9, 3.141592, 2.718]
    XCTAssertEqual(operandA.count, 7)
    XCTAssertEqual(operandB.count, 7)
    
    // Encode the GPU command.
    for i in operandA.indices {
      contents[i] = operandA[i]
      contents[8 + i] = operandB[i]
    }
    let encoder = commandBuffer.makeComputeCommandEncoder()!
    encoder.setComputePipelineState(pipeline)
    encoder.setBuffer(buffer, offset: 0, index: 0)
    encoder.setBuffer(buffer, offset: 32, index: 1)
    encoder.setBuffer(buffer, offset: 64, index: 2)
    encoder.dispatchThreads(
      MTLSizeMake(7, 1, 1), threadsPerThreadgroup: MTLSizeMake(128, 1, 1))
    encoder.endEncoding()
    commandBuffer.commit()
    commandBuffer.waitUntilCompleted()
    
    // Compute the expected result.
    var operandC: [Float] = []
    for i in operandA.indices {
      let valueA = operandA[i]
      let valueB = operandB[i]
      let valueC = valueA + valueB
      operandC.append(valueC)
    }
    
    // Check that the CPU and GPU results match to machine precision.
    for i in operandA.indices {
      let expected = operandC[i]
      let actual = contents[16 + i]
      XCTAssertEqual(expected, actual)
    }
  }
  
  // Test that async copies compile and execute properly at runtime.
  func testAsyncCopy() throws {
    
  }
}
