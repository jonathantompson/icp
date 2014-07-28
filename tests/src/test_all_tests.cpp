//
//  test_all_tests.cpp
//
//  Created by Jonathan Tompson on 2/23/13.
//
//  Test everything

//#include "test_callback.h"
//#include "test_callback_queue.h"
//#include "test_data_str.h"
//#include "test_math.h"
//#include "test_thread.h"
//#include "test_thread_pool.h"
//#include "test_optimization.h"
#include "test_iterative_closest_point.h"
//#include "test_math/test_profile_simd_math.h"  // Profile last

#include "icp/debug_util/debug_util.h"  // Must come last in .cpp with main

using std::cout;
using std::endl;

int main(int argc, char* argv[]) {
#if defined(_DEBUG) || defined(DEBUG)
  icp::debug::EnableMemoryLeakChecks();
  // icp::debug::SetBreakPointOnAlocation(26730);
#endif
  int ret_val = RUN_TESTS(argc, argv);
#ifdef _WIN32
  system("PAUSE");
#endif
  return ret_val;
}
