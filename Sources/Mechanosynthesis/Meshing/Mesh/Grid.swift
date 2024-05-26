//
//  Grid.swift
//
//
//  Created by Philip Turner on 5/26/24.
//

public struct Grid<Element> {
  // The lower corner of the bounding box enclosing the cells.
  var offset: SIMD3<UInt32>
  
  // The size of the bounding box enclosing the cells.
  var dimensions: SIMD3<UInt32>
  
  // A uniform grid of elements.
  var cells: [Element]
}
