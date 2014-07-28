#include <string>
#ifdef __APPLE__
  #include <signal.h>
#endif
#include "icp/exceptions/wruntime_error.h"
#include "icp/string_util/string_util.h"

using icp::string_util::ToNarrowString;
using icp::string_util::ToWideString;

#if (defined(DEBUG) || defined(_DEBUG)) && defined(BREAK_ON_EXCEPTION)
  #define BREAK_ON_EXCEPTION_INT
#endif

namespace std {

  wruntime_error::wruntime_error(const wstring& errorMsg)
    : runtime_error(ToNarrowString(errorMsg)), mErrorMsg(errorMsg) {
    // NOTE: We give the runtime_error base the narrow version of the 
    //  error message. This is what will get shown if what() is called.
    //  The wruntime_error inserter or errorMsg() should be used to get 
    //  the wide version.
#if defined(_WIN32) && defined(BREAK_ON_EXCEPTION_INT)
    if (IsDebuggerPresent()) {
      _CrtDbgBreak();
    }
#endif
#if defined(__APPLE__) && defined(BREAK_ON_EXCEPTION_INT)
    raise(SIGTRAP);
    //__builtin_trap();
#endif
  }

  wruntime_error::wruntime_error(const string& errorMsg)
    : runtime_error(errorMsg), mErrorMsg(ToWideString(errorMsg)) {
    // NOTE: We give the runtime_error base the narrow version of the 
    //  error message. This is what will get shown if what() is called.
    //  The wruntime_error inserter or errorMsg() should be used to get 
    //  the wide version.
#if defined(_WIN32) && defined(BREAK_ON_EXCEPTION_INT)
    if (IsDebuggerPresent()) {
      _CrtDbgBreak();
    }
#endif
#if defined(__APPLE__) && defined(BREAK_ON_EXCEPTION_INT)
      raise(SIGTRAP);
      //__builtin_trap();
#endif
  }
   
  wruntime_error::wruntime_error(const wruntime_error& rhs)
    : runtime_error(ToNarrowString(rhs.errorMsg())), mErrorMsg(rhs.errorMsg()) {
  }

  wruntime_error& wruntime_error::operator=(const wruntime_error& rhs) {
    // copy the wruntime_error
    runtime_error::operator=(rhs); 
    mErrorMsg = rhs.mErrorMsg; 

#if defined(_WIN32) && defined(BREAK_ON_EXCEPTION_INT)
    if (IsDebuggerPresent()) {
      _CrtDbgBreak();
    }
#endif
#if defined(__APPLE__) && defined(BREAK_ON_EXCEPTION_INT)
    raise(SIGTRAP);
    //__builtin_trap();
#endif
    return *this; 
  }

#if defined(WIN32) || defined(_WIN32) || defined(__APPLE__)
  wruntime_error::~wruntime_error() {
  }
#else
  wruntime_error::~wruntime_error() _GLIBCXX_USE_NOEXCEPT {
  }
#endif

  const wstring& wruntime_error::errorMsg() const { 
    return mErrorMsg; 
  }
}
