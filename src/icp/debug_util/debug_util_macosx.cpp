#ifdef __APPLE__

#include <stdio.h>  // printf
#include "icp/debug_util/debug_util.h"

namespace icp {
namespace debug {

  void EnableMemoryLeakChecks() {
    printf("WARNING: EnableMemoryLeakChecks not implemented yet for Mac\n");
  }

  void SetBreakPointOnAlocation(int alloc_num) {
    static_cast<void>(alloc_num);
    printf("WARNING: SetBreakPointOnAlocation not implemented yet for Mac\n");
  }

}  // namespace debug
}  // namespace icp

#endif  // __APPLE__
