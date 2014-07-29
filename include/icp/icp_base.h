//
//  icp.h
//
//  Created by Jonathan Tompson on 2/23/13.
//
//  Top level header file.  Includes most of the useful utilities.
//  Everything is in the icp namespace!
//  

#pragma once

#include "icp/math/math_types.h"  // Lots of default matrix and vector types
#include "icp/math/noise.h"  // Generate continuously varying noisy keyframes
#include "icp/string_util/string_util.h"  // Common string functions
#include "icp/exceptions/wruntime_error.h"  // Breakpoint in DEBUG builds
#if defined( _WIN32 )
  #include "icp/string_util/win32_debug_buffer.h"
#endif
#include "icp/math/perlin_noise.h"
#include "icp/math/icp.h"
