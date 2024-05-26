//
//  Grid.swift
//
//
//  Created by Philip Turner on 5/26/24.
//

public struct GridDescriptor<Element> {
  // The lower corner of the bounding box enclosing the cells.
  public var offset: SIMD3<UInt32> = .zero
  
  // The size of the bounding box enclosing the cells.
  public var dimensions: SIMD3<UInt32>?
  
  // An empty element to fill the grid with.
  public var emptyElement: Element?
  
  public init() {
    
  }
}

public struct Grid<Element> {
  // The lower corner of the bounding box enclosing the cells.
  public var offset: SIMD3<UInt32> = .zero
  
  // The size of the bounding box enclosing the cells.
  public var dimensions: SIMD3<UInt32>?
  
  // A uniform grid of elements.
  public var cells: [Element]
  
  public init(descriptor: GridDescriptor<Element>) {
    guard let dimensions = descriptor.dimensions,
          let emptyElement = descriptor.emptyElement else {
      fatalError("Descriptor was incomplete.")
    }
    guard all(dimensions .>= 0) else {
      fatalError("Cannot have negative dimensions.")
    }
    
    // Allocate an array for the cells.
    let cellCount = dimensions[0] * dimensions[1] * dimensions[2]
    cells = Array(repeating: emptyElement, count: Int(cellCount))
  }
}
