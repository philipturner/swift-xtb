import XCTest
import Mechanosynthesis
import Numerics

final class OctreeTests: XCTestCase {
  func testOctree() throws {
    var descriptor = OctreeDescriptor()
    descriptor.sizeExponent = 4
    var octree = Octree(descriptor: descriptor)
    XCTAssertEqual(octree.nodes.count, 1)
    XCTAssertEqual(octree.nodes[0], OctreeNode(
      centerAndSpacing: SIMD4(0, 0, 0, 16),
      parentIndex: UInt32.max,
      branchesIndex: UInt32.max,
      branchesMask: SIMD8.zero))
    
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
      branchesMask: SIMD8(0, 1, 0, 0, 0, 0, 0, 0)))
    XCTAssertEqual(octree.nodes[1], OctreeNode(
      centerAndSpacing: SIMD4(4, -4, -4, 8),
      parentIndex: 0,
      branchesIndex: UInt32.max,
      branchesMask: SIMD8.zero))
    
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
      branchesMask: SIMD8(1, 1, 0, 0, 0, 0, 0, 0)))
    XCTAssertEqual(octree.nodes[1], OctreeNode(
      centerAndSpacing: SIMD4(-4, -4, -4, 8),
      parentIndex: 0,
      branchesIndex: UInt32.max,
      branchesMask: SIMD8.zero))
    XCTAssertEqual(octree.nodes[2], OctreeNode(
      centerAndSpacing: SIMD4(4, -4, -4, 8),
      parentIndex: 0,
      branchesIndex: 3,
      branchesMask: SIMD8(0, 0, 0, 0, 1, 0, 0, 0)))
    XCTAssertEqual(octree.nodes[3], OctreeNode(
      centerAndSpacing: SIMD4(2, -6, -2, 4),
      parentIndex: 2,
      branchesIndex: UInt32.max,
      branchesMask: SIMD8.zero))
    
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
      branchesMask: SIMD8(1, 1, 0, 0, 0, 0, 0, 0)))
    XCTAssertEqual(octree.nodes[1], OctreeNode(
      centerAndSpacing: SIMD4(-4, -4, -4, 8),
      parentIndex: 0,
      branchesIndex: UInt32.max,
      branchesMask: SIMD8.zero))
    XCTAssertEqual(octree.nodes[2], OctreeNode(
      centerAndSpacing: SIMD4(4, -4, -4, 8),
      parentIndex: 0,
      branchesIndex: UInt32.max,
      branchesMask: SIMD8.zero))
    
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
      branchesMask: SIMD8(0, 1, 0, 0, 0, 0, 0, 0)))
    XCTAssertEqual(octree.nodes[1], OctreeNode(
      centerAndSpacing: SIMD4(4, -4, -4, 8),
      parentIndex: 0,
      branchesIndex: 2,
      branchesMask: SIMD8(0, 0, 0, 0, 1, 0, 0, 0)))
    XCTAssertEqual(octree.nodes[2], OctreeNode(
      centerAndSpacing: SIMD4(2, -6, -2, 4),
      parentIndex: 1,
      branchesIndex: UInt32.max,
      branchesMask: SIMD8.zero))
  }
}
