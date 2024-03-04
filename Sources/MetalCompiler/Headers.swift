//
//  Headers.swift
//
//
//  Created by Philip Turner on 3/3/24.
//

let metalSimdgroupEvent: String = """
// -*- Metal -*-
//===-- metal_simdgroup_event ---------------------------------------------===//
// Copyright (c) 2023 Philip Turner. See MIT LICENSE
//===----------------------------------------------------------------------===//

#ifndef __METAL_SIMDGROUP_EVENT
#define __METAL_SIMDGROUP_EVENT

// Contains C++ symbols accessible to a developer through automatic code
// completion in Xcode 14.2. Formatted with the same style as the Metal Standard
// Library for consistency with other Metal code.

#if defined(__HAVE_SIMDGROUP_FUTURE__)
#include <metal_simdgroup_future>
#endif

#pragma METAL internals : enable
namespace metal
{
#if !defined(__HAVE_SIMDGROUP_FUTURE__)
  enum class simdgroup_async_copy_clamp_mode {
    clamp_to_zero = 0,
    clamp_to_edge = 1
  };
#endif
  
  struct simdgroup_event {
    METAL_FUNC simdgroup_event() thread {}
    
    template <typename T>
    METAL_FUNC void async_copy(threadgroup T *dst, const device T *src, ulong n_elements) thread {
#if defined(__HAVE_SIMDGROUP_FUTURE__)
      event = __metal_simdgroup_async_copy_1d(sizeof(T), alignof(T), reinterpret_cast<threadgroup void *>(dst), reinterpret_cast<const device void *>(src), n_elements);
#endif
    }
    
    template <typename T>
    METAL_FUNC void async_copy(device T *dst, const threadgroup T *src, ulong n_elements) thread {
#if defined(__HAVE_SIMDGROUP_FUTURE__)
      event = __metal_simdgroup_async_copy_1d(sizeof(T), alignof(T), reinterpret_cast<device void *>(dst), reinterpret_cast<const threadgroup void *>(src), n_elements);
#endif
    }
    
    template <typename T>
    METAL_FUNC void async_copy(threadgroup T *dst, ushort dst_elements_per_row, ushort2 dst_tile_dimensions, const device T *src, uint src_elements_per_row, ushort2 src_tile_dimensions, bool transpose_matrix = false, simdgroup_async_copy_clamp_mode clamp_mode = simdgroup_async_copy_clamp_mode::clamp_to_zero) thread {
#if defined(__HAVE_SIMDGROUP_FUTURE__)
      if (transpose_matrix) {
        src_tile_dimensions = src_tile_dimensions.yx;
        dst_tile_dimensions = dst_tile_dimensions.yx;
      }
      event = __metal_simdgroup_async_copy_2d(sizeof(T), alignof(T), reinterpret_cast<threadgroup void *>(dst), ushort(dst_elements_per_row), 1, ulong2(dst_tile_dimensions), reinterpret_cast<const device void *>(src), uint(src_elements_per_row), 1, ulong2(src_tile_dimensions), long2(0), static_cast<int>(clamp_mode));
#endif
    }
    
    template <typename T>
    METAL_FUNC void async_copy(device T *dst, uint dst_elements_per_row, ushort2 dst_tile_dimensions, const threadgroup T *src, ushort src_elements_per_row, ushort2 src_tile_dimensions, bool transpose_matrix = false) thread {
#if defined(__HAVE_SIMDGROUP_FUTURE__)
      if (transpose_matrix) {
        src_tile_dimensions = src_tile_dimensions.yx;
        dst_tile_dimensions = dst_tile_dimensions.yx;
      }
      event = __metal_simdgroup_async_copy_2d(sizeof(T), alignof(T), reinterpret_cast<device void *>(dst), uint(dst_elements_per_row), 1, ulong2(dst_tile_dimensions), reinterpret_cast<const threadgroup void *>(src), ushort(src_elements_per_row), 1, ulong2(src_tile_dimensions), long2(0), 0);
#endif
    }
    
    METAL_FUNC static void wait(int count, thread simdgroup_event *events) {
#if defined(__HAVE_SIMDGROUP_FUTURE__)
      __metal_wait_simdgroup_events(count, reinterpret_cast<thread __metal_simdgroup_event_t*>(events));
#endif
    }
    
  private:
#if defined(__HAVE_SIMDGROUP_FUTURE__)
    __metal_simdgroup_event_t event;
#endif
  };
} // namespace metal
#pragma METAL internals : disable

#endif // __METAL_SIMDGROUP_EVENT
"""

let metalSimdgroupMatrixStorage: String = """
// -*- Metal -*-
//===-- metal_simdgroup_matrix_storage ------------------------------------===//
// Copyright (c) 2023 Philip Turner. See MIT LICENSE
//===----------------------------------------------------------------------===//

#ifndef __METAL_SIMDGROUP_MATRIX_STORAGE
#define __METAL_SIMDGROUP_MATRIX_STORAGE

// Contains C++ symbols accessible to a developer through automatic code
// completion in Xcode 14.2. Formatted with the same style as the Metal Standard
// Library for consistency with other Metal code.

#if defined(__HAVE_SIMDGROUP_MATRIX__)
#pragma METAL internals : enable
namespace metal
{
  template <typename T>
  struct simdgroup_matrix_storage {
    typedef vec<T, 64> storage_type;
    
    storage_type t;
    
    METAL_FUNC thread vec<T, 2>* thread_elements() thread {
      return reinterpret_cast<thread vec<T, 2>*>(&t);
    }
    
    METAL_FUNC simdgroup_matrix_storage() thread = default;
    
    METAL_FUNC simdgroup_matrix_storage(vec<T, 2> thread_elements) thread {
      *(this->thread_elements()) = thread_elements;
    }
    
    METAL_FUNC static ushort2 offset(ushort thread_index_in_simdgroup) {
      // https://patents.google.com/patent/US11256518B2
      ushort lane_id = thread_index_in_simdgroup;
      ushort quad_id = lane_id / 4;
      
      constexpr ushort QUADRANT_SPAN_M = 4;
      constexpr ushort THREADS_PER_QUADRANT = 8;
      ushort M_floor_of_quadrant = (quad_id / 4) * QUADRANT_SPAN_M;
      ushort M_in_quadrant = (lane_id / 2) % (THREADS_PER_QUADRANT / 2);
      ushort M_in_simd = M_floor_of_quadrant + M_in_quadrant;
      
      ushort N_floor_of_quadrant = (quad_id & 2) * 2; // 0 or 4
      ushort N_in_quadrant = (lane_id % 2) * 2; // 0 or 2
      ushort N_in_simd = N_floor_of_quadrant + N_in_quadrant;
      
      return ushort2(N_in_simd, M_in_simd);
    }
    
    METAL_FUNC static device T* apply_offset(device T *src, uint elements_per_row, uint2 matrix_origin, bool transpose_matrix = false) {
      if (transpose_matrix) {
        return src + ulong(matrix_origin.x * elements_per_row) + matrix_origin.y;
      } else {
        return src + ulong(matrix_origin.y * elements_per_row) + matrix_origin.x;
      }
    }
    
    METAL_FUNC static threadgroup T* apply_offset(threadgroup T *src, ushort elements_per_row, ushort2 matrix_origin, bool transpose_matrix = false) {
      if (transpose_matrix) {
        return src + matrix_origin.x * elements_per_row + matrix_origin.y;
      } else {
        return src + matrix_origin.y * elements_per_row + matrix_origin.x;
      }
    }
    
    // WARNING: All load and store functions assume the X dimension is divisible by 2.
    
    METAL_FUNC void load(const device T *src, uint elements_per_row, ushort2 matrix_origin, bool transpose_matrix = false) {
      if (transpose_matrix) {
        *(thread_elements()) = vec<T, 2>(src[ulong(matrix_origin.x * elements_per_row) + matrix_origin.y], src[ulong((matrix_origin.x + 1) * elements_per_row) + matrix_origin.y]);
      } else {
        *(thread_elements()) = *reinterpret_cast<const device vec<T, 2>*>(src + ulong(matrix_origin.y * elements_per_row) + matrix_origin.x);
      }
    }
    
    METAL_FUNC void load(const threadgroup T *src, ushort elements_per_row, ushort2 matrix_origin, bool transpose_matrix = false) {
      if (transpose_matrix) {
        *(thread_elements()) = vec<T, 2>(src[matrix_origin.x * elements_per_row + matrix_origin.y], src[(matrix_origin.x + 1) * elements_per_row + matrix_origin.y]);
      } else {
        *(thread_elements()) = *reinterpret_cast<const threadgroup vec<T, 2>*>(src + matrix_origin.y * elements_per_row + matrix_origin.x);
      }
    }
    
    METAL_FUNC void load_first(const device T *src, ushort elements_per_row, ushort2 matrix_origin, bool transpose_matrix = false) {
      if (transpose_matrix) {
        thread_elements()[0][0] = src[matrix_origin.x * elements_per_row + matrix_origin.y];
      } else {
        thread_elements()[0][0] = src[matrix_origin.y * elements_per_row + matrix_origin.x];
      }
    }
    
    METAL_FUNC void load_second(const device T *src, ushort elements_per_row, ushort2 matrix_origin, bool transpose_matrix = false) {
      if (transpose_matrix) {
        thread_elements()[0][1] = src[matrix_origin.x * elements_per_row + matrix_origin.y];
      } else {
        thread_elements()[0][1] = src[matrix_origin.y * elements_per_row + matrix_origin.x];
      }
    }
    
    METAL_FUNC void store(device T *dst, uint elements_per_row, ushort2 matrix_origin, bool transpose_matrix = false) {
      if (transpose_matrix) {
        dst[ulong(matrix_origin.x * elements_per_row) + matrix_origin.y] = thread_elements()[0][0];
        dst[ulong((matrix_origin.x + 1) * elements_per_row) + matrix_origin.y] = thread_elements()[0][1];
      } else {
        *reinterpret_cast<device vec<T, 2>*>(dst + matrix_origin.y * elements_per_row + matrix_origin.x) = *(thread_elements());
      }
    }
    
    METAL_FUNC void store_first(device T *dst, uint elements_per_row, ushort2 matrix_origin, bool transpose_matrix = false) {
      if (transpose_matrix) {
        dst[ulong(matrix_origin.x * elements_per_row) + matrix_origin.y] = thread_elements()[0][0];
      } else {
        dst[matrix_origin.y * elements_per_row + matrix_origin.x] = thread_elements()[0][0];
      }
    }
    
    METAL_FUNC void store_second(device T *dst, uint elements_per_row, ushort2 matrix_origin, bool transpose_matrix = false) {
      if (transpose_matrix) {
        dst[ulong(matrix_origin.x * elements_per_row) + matrix_origin.y] = thread_elements()[0][1];
      } else {
        dst[matrix_origin.y * elements_per_row + matrix_origin.x] = thread_elements()[0][1];
      }
    }
    
    METAL_FUNC void store(threadgroup T *dst, ushort elements_per_row, ushort2 matrix_origin, bool transpose_matrix = false) {
      if (transpose_matrix) {
        dst[matrix_origin.x * elements_per_row + matrix_origin.y] = thread_elements()[0][0];
        dst[(matrix_origin.x + 1) * elements_per_row + matrix_origin.y] = thread_elements()[0][1];
      } else {
        *reinterpret_cast<threadgroup vec<T, 2>*>(dst + matrix_origin.y * elements_per_row + matrix_origin.x) = *(thread_elements());
      }
    }
    
    template <typename U, typename V>
    METAL_FUNC void multiply(simdgroup_matrix_storage<U> a, simdgroup_matrix_storage<V> b, bool accumulate = true) {
      if (!accumulate) {
        *(thread_elements()) = vec<T, 2>(0);
      }
      t = __metal_simdgroup_matrix_8x8_multiply_accumulate(a.t, b.t, t, typename simdgroup_matrix_storage<T>::storage_type());
    }
  };
} // namespace metal
#pragma METAL internals : disable
#endif

#endif // __METAL_SIMDGROUP_MATRIX_STORAGE
"""
