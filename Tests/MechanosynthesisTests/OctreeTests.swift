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
      branchesMask: 0b0000_0000))
    
    octree.resizeNodes(expanded: [(0, SIMD8(0, 0b0000_0010, 0, 0, 0, 0, 0, 0))], contracted: [])
    XCTAssertEqual(octree.nodes.count, 2)
    XCTAssertEqual(octree.nodes[0], OctreeNode(
      centerAndSpacing: SIMD4(0, 0, 0, 16),
      parentIndex: UInt32.max,
      branchesIndex: 1,
      branchesMask: 0b0000_0010))
    XCTAssertEqual(octree.nodes[1], OctreeNode(
      centerAndSpacing: SIMD4(4, -4, -4, 8),
      parentIndex: 0,
      branchesIndex: UInt32.max,
      branchesMask: 0b0000_0000))
    
    octree.resizeNodes(
      expanded: [(0, SIMD8(0b0000_0001, 0, 0, 0, 0, 0, 0, 0)), (1, SIMD8(0, 0, 0, 0, 0b0001_0000, 0, 0, 0))], contracted: [])
    XCTAssertEqual(octree.nodes.count, 4)
    XCTAssertEqual(octree.nodes[0], OctreeNode(
      centerAndSpacing: SIMD4(0, 0, 0, 16),
      parentIndex: UInt32.max,
      branchesIndex: 1,
      branchesMask: 0b0000_0011))
    XCTAssertEqual(octree.nodes[1], OctreeNode(
      centerAndSpacing: SIMD4(-4, -4, -4, 8),
      parentIndex: 0,
      branchesIndex: UInt32.max,
      branchesMask: 0b0000_0000))
    XCTAssertEqual(octree.nodes[2], OctreeNode(
      centerAndSpacing: SIMD4(4, -4, -4, 8),
      parentIndex: 0,
      branchesIndex: 3,
      branchesMask: 0b0001_0000))
    XCTAssertEqual(octree.nodes[3], OctreeNode(
      centerAndSpacing: SIMD4(2, -6, -2, 4),
      parentIndex: 2,
      branchesIndex: UInt32.max,
      branchesMask: 0b0000_0000))
    
    octree.resizeNodes(expanded: [], contracted: [3])
    XCTAssertEqual(octree.nodes.count, 3)
    XCTAssertEqual(octree.nodes[0], OctreeNode(
      centerAndSpacing: SIMD4(0, 0, 0, 16),
      parentIndex: UInt32.max,
      branchesIndex: 1,
      branchesMask: 0b0000_0011))
    XCTAssertEqual(octree.nodes[1], OctreeNode(
      centerAndSpacing: SIMD4(-4, -4, -4, 8),
      parentIndex: 0,
      branchesIndex: UInt32.max,
      branchesMask: 0b0000_0000))
    XCTAssertEqual(octree.nodes[2], OctreeNode(
      centerAndSpacing: SIMD4(4, -4, -4, 8),
      parentIndex: 0,
      branchesIndex: UInt32.max,
      branchesMask: 0b0000_0000))
    
    octree.resizeNodes(
      expanded: [(2, SIMD8(0, 0, 0, 0, 0b0001_0000, 0, 0, 0))], contracted: [1])
    XCTAssertEqual(octree.nodes.count, 3)
    XCTAssertEqual(octree.nodes[0], OctreeNode(
      centerAndSpacing: SIMD4(0, 0, 0, 16),
      parentIndex: UInt32.max,
      branchesIndex: 1,
      branchesMask: 0b0000_0010))
    XCTAssertEqual(octree.nodes[1], OctreeNode(
      centerAndSpacing: SIMD4(4, -4, -4, 8),
      parentIndex: 0,
      branchesIndex: 2,
      branchesMask: 0b0001_0000))
    XCTAssertEqual(octree.nodes[2], OctreeNode(
      centerAndSpacing: SIMD4(2, -6, -2, 4),
      parentIndex: 1,
      branchesIndex: UInt32.max,
      branchesMask: 0b0000_0000))
  }
}
