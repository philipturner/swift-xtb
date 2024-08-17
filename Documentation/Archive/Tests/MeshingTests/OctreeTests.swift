import XCTest
import Meshing

final class OctreeTests: XCTestCase {
  func testPrefixSum() throws {
    func prefixSum(_ x: SIMD8<UInt8>) -> SIMD8<UInt8> {
      // 00000000
      // 11111111
      var output: UInt64 = .zero
      let forwardSum0 = unsafeBitCast(x, to: UInt64.self)
      
      // 01010101
      // 02020202
      let shifted1 = (forwardSum0 & 0x00FF00FF00FF00FF) << 8
      let forwardSum1 = (forwardSum0 &+ shifted1) & 0xFF00FF00FF00FF00
      output = shifted1
      
      // 01230123
      // 00040004
      let shifted2 = (forwardSum1 & 0x0000FF000000FF00) << 16
      let forwardSum2 = (forwardSum1 &+ shifted2) & 0xFF000000FF000000
      output &+= (forwardSum1 & 0x0000FF000000FF00) &* (0x0101 << 8)
      
      // 01234567
      output &+= forwardSum2 &* (0x01010101) << 8
      return unsafeBitCast(output, to: SIMD8<UInt8>.self)
    }
    
    do {
      var x = SIMD8<UInt8>(1, 1, 1, 1, 1, 1, 1, 1)
      x = prefixSum(x)
      XCTAssertEqual(x, SIMD8(0, 1, 2, 3, 4, 5, 6, 7))
    }
    
    do {
      var x = SIMD8<UInt8>(1, 9, 1, 2, 1, 1, 1, 1)
      x = prefixSum(x)
      XCTAssertEqual(x, SIMD8(0, 1, 10, 11, 13, 14, 15, 16))
    }
    
    do {
      let branchesMask = SIMD8<UInt8>(0, 0, 255, 255, 255, 255, 255, 255)
      var x = SIMD8<UInt8>(repeating: 0)
      x.replace(with: 1, where: branchesMask .!= 255)
      x = prefixSum(x)
      x.replace(with: 255, where: branchesMask .== 255)
      XCTAssertEqual(x, SIMD8(0, 1, 255, 255, 255, 255, 255, 255))
    }
    
    do {
      let branchesMask = SIMD8<UInt8>(255, 0, 255, 255, 0, 0, 255, 255)
      var x = SIMD8<UInt8>(repeating: 0)
      x.replace(with: 1, where: branchesMask .!= 255)
      x = prefixSum(x)
      x.replace(with: 255, where: branchesMask .== 255)
      XCTAssertEqual(x, SIMD8(255, 0, 255, 255, 1, 2, 255, 255))
    }
  }
  
  func testOctree() throws {
    var descriptor = OctreeDescriptor()
    descriptor.sizeExponent = 4
    var octree = Octree(descriptor: descriptor)
    XCTAssertEqual(octree.nodes.count, 1)
    XCTAssertEqual(octree.nodes[0], OctreeNode(
      centerAndSpacing: SIMD4(0, 0, 0, 16),
      parentIndex: UInt32.max,
      branchesIndex: UInt32.max,
      branchesMask: SIMD8(repeating: UInt8.max)))
    
    octree.resizeNodes(
      expanded: [
        SIMD8(0, 1, 0, 0, 0, 0, 0, 0),
      ],
      contracted: [
        false,
      ])
    XCTAssertEqual(octree.nodes.count, 2)
    XCTAssertEqual(octree.nodes[0], OctreeNode(
      centerAndSpacing: SIMD4(0, 0, 0, 16),
      parentIndex: UInt32.max,
      branchesIndex: 1,
      branchesMask: SIMD8(255, 0, 255, 255, 255, 255, 255, 255)))
    XCTAssertEqual(octree.nodes[1], OctreeNode(
      centerAndSpacing: SIMD4(4, -4, -4, 8),
      parentIndex: 0,
      branchesIndex: UInt32.max,
      branchesMask: SIMD8(repeating: UInt8.max)))
    
    octree.resizeNodes(
      expanded: [
        SIMD8(1, 0, 0, 0, 0, 0, 0, 0),
        SIMD8(0, 0, 0, 0, 1, 0, 0, 0),
      ],
      contracted: [
        false,
        false,
      ])
    XCTAssertEqual(octree.nodes.count, 4)
    XCTAssertEqual(octree.nodes[0], OctreeNode(
      centerAndSpacing: SIMD4(0, 0, 0, 16),
      parentIndex: UInt32.max,
      branchesIndex: 1,
      branchesMask: SIMD8(0, 1, 255, 255, 255, 255, 255, 255)))
    XCTAssertEqual(octree.nodes[1], OctreeNode(
      centerAndSpacing: SIMD4(-4, -4, -4, 8),
      parentIndex: 0,
      branchesIndex: UInt32.max,
      branchesMask: SIMD8(repeating: UInt8.max)))
    XCTAssertEqual(octree.nodes[2], OctreeNode(
      centerAndSpacing: SIMD4(4, -4, -4, 8),
      parentIndex: 0,
      branchesIndex: 3,
      branchesMask: SIMD8(255, 255, 255, 255, 0, 255, 255, 255)))
    XCTAssertEqual(octree.nodes[3], OctreeNode(
      centerAndSpacing: SIMD4(2, -6, -2, 4),
      parentIndex: 2,
      branchesIndex: UInt32.max,
      branchesMask: SIMD8(repeating: UInt8.max)))
    
    octree.resizeNodes(
      expanded: [
        SIMD8.zero,
        SIMD8.zero,
        SIMD8.zero,
      ],
      contracted: [
        false,
        false,
        false,
        true,
      ])
    XCTAssertEqual(octree.nodes.count, 3)
    XCTAssertEqual(octree.nodes[0], OctreeNode(
      centerAndSpacing: SIMD4(0, 0, 0, 16),
      parentIndex: UInt32.max,
      branchesIndex: 1,
      branchesMask: SIMD8(0, 1, 255, 255, 255, 255, 255, 255)))
    XCTAssertEqual(octree.nodes[1], OctreeNode(
      centerAndSpacing: SIMD4(-4, -4, -4, 8),
      parentIndex: 0,
      branchesIndex: UInt32.max,
      branchesMask: SIMD8(repeating: UInt8.max)))
    XCTAssertEqual(octree.nodes[2], OctreeNode(
      centerAndSpacing: SIMD4(4, -4, -4, 8),
      parentIndex: 0,
      branchesIndex: UInt32.max,
      branchesMask: SIMD8(repeating: UInt8.max)))
    
    octree.resizeNodes(
      expanded: [
        SIMD8.zero,
        SIMD8.zero,
        SIMD8(0, 0, 0, 0, 1, 0, 0, 0),
        SIMD8.zero,
      ],
      contracted: [
        false,
        true,
        false,
      ])
    XCTAssertEqual(octree.nodes.count, 3)
    XCTAssertEqual(octree.nodes[0], OctreeNode(
      centerAndSpacing: SIMD4(0, 0, 0, 16),
      parentIndex: UInt32.max,
      branchesIndex: 1,
      branchesMask: SIMD8(255, 0, 255, 255, 255, 255, 255, 255)))
    XCTAssertEqual(octree.nodes[1], OctreeNode(
      centerAndSpacing: SIMD4(4, -4, -4, 8),
      parentIndex: 0,
      branchesIndex: 2,
      branchesMask: SIMD8(255, 255, 255, 255, 0, 255, 255, 255)))
    XCTAssertEqual(octree.nodes[2], OctreeNode(
      centerAndSpacing: SIMD4(2, -6, -2, 4),
      parentIndex: 1,
      branchesIndex: UInt32.max,
      branchesMask: SIMD8(repeating: UInt8.max)))
  }
}
