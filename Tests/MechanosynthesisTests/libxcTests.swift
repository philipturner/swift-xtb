import XCTest
import libxc

// Reproduction of tests from the LibXC 5.1.X manual:
// https://www.tddft.org/programs/libxc/manual/libxc-5.1.x
final class libxcTests: XCTestCase {
  func testSmallProgram() throws {
    var rho: [Double] = [0.1, 0.2, 0.3, 0.4, 0.5]
    var sigma: [Double] = [0.2, 0.3, 0.4, 0.5, 0.6]
    var exc = [Double](repeating: -1, count: 5)
    
    // Get the LibXC version.
    var vmajor: Int32 = -1
    var vminor: Int32 = -1
    var vmicro: Int32 = -1
    xc_version(&vmajor, &vminor, &vmicro)
    XCTAssertGreaterThanOrEqual(vmajor, 6)
    XCTAssertGreaterThanOrEqual(vminor, 0)
    XCTAssertGreaterThanOrEqual(vmicro, 0)
    
    // Initialize the functional.
    var functional = xc_func_type()
    let error = xc_func_init(&functional, XC_LDA_X, XC_UNPOLARIZED)
    XCTAssertEqual(error, 0)
    
    // Evaluate the energy density, depending on the family.
    switch functional.info.pointee.family {
    case XC_FAMILY_LDA:
      xc_lda_exc(&functional, 5, &rho, &exc)
    case XC_FAMILY_GGA:
      fallthrough
    case XC_FAMILY_HYB_GGA:
      xc_gga_exc(&functional, 5, &rho, &sigma, &exc)
    default:
      fatalError("Unexpected family: \(functional.info.pointee.family)")
    }
    
    // Print out density and energy density per particle.
    XCTAssertEqual(exc[0], -0.342809, accuracy: 1e-6)
    XCTAssertEqual(exc[1], -0.431912, accuracy: 1e-6)
    XCTAssertEqual(exc[2], -0.494416, accuracy: 1e-6)
    XCTAssertEqual(exc[3], -0.544175, accuracy: 1e-6)
    XCTAssertEqual(exc[4], -0.586194, accuracy: 1e-6)
    
    // Deallocate memory.
    xc_func_end(&functional)
  }
  
  func testSlaterExchange() throws {
    let name: StaticString = "Slater exchange"
    let flags = XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC |
    XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC
    
    let xc_func_info_lda_x = xc_func_info_type(
      number: XC_LDA_X,
      kind: XC_EXCHANGE,
      name: UnsafePointer<Int8>(OpaquePointer(name.utf8Start)),
      family: XC_FAMILY_LDA,
      
      // (xc_ref_Dirac1930_376, xc_ref_Bloch_1929_545, nil, nil, nil)
      refs: (nil, nil, nil, nil, nil),
      flags: flags,
      dens_threshold: 1e-24,
      ext_params: func_params_type(),
      init: nil,
      end: nil,
      lda: nil,
      gga: nil,
      mgga: nil)
    _ = xc_func_info_lda_x
  }
  
  func testInfoStructure() throws {
    var output: String = ""
    
    // Initialize the functional.
    var functional = xc_func_type()
    xc_func_init(&functional, XC_GGA_X_B88, XC_UNPOLARIZED)
    
    // Print out information.
    output += "The functional '"
    output += String(cString: functional.info.pointee.name)
    output += "' is "
    switch functional.info.pointee.kind {
    case XC_EXCHANGE:
      output += "an exchange functional"
    case XC_CORRELATION:
      output += "a correlation functional"
    case XC_EXCHANGE_CORRELATION:
      output += "an exchange-correlation functional"
    case XC_KINETIC:
      output += "a kinetic energy functional"
    default:
      output += "of unknown kind"
    }
    output += ", it belongs to the '"
    
    // Print out family.
    switch functional.info.pointee.family {
    case XC_FAMILY_LDA:
      output += "LDA"
    case XC_FAMILY_GGA:
      output += "GGA"
    case XC_FAMILY_HYB_GGA:
      output += "Hybrid GGA"
    case XC_FAMILY_MGGA:
      output += "MGGA"
    case XC_FAMILY_HYB_MGGA:
      output += "Hybrid MGGA"
    default:
      output += "unknown"
    }
    output += "' family and is defined in the reference(s):\n"
    
    // Print out references.
    let refs = functional.info.pointee.refs
    let refsArray = [refs.0, refs.1, refs.2, refs.3, refs.4]
    var referenceID = 0
    for ref in refsArray where ref != nil {
      referenceID += 1
      output += "[\(referenceID)] "
      output += String(cString: ref!.pointee.ref)
      output += "\n"
    }
    
    // Deallocate memory.
    xc_func_end(&functional)
    
    XCTAssertEqual(output, """
      The functional 'Becke 88' is an exchange functional, it belongs to the \
      'GGA' family and is defined in the reference(s):
      [1] A. D. Becke.,  Phys. Rev. A 38, 3098 (1988)
      
      """)
  }
}
