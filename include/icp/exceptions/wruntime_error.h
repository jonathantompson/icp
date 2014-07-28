//
//  wruntime_error.h
//
//  Created by Jonathan Tompson on 5/27/12.
//
//  A simple wide version of std::runtime_error exception with break points for
//  visual studio.

#pragma once

#include <stdexcept>  // for runtime_error
#include <string>

namespace std {

  class wruntime_error : public std::runtime_error {
  public:       
    explicit wruntime_error(const std::string& errorMsg);
    explicit wruntime_error(const std::wstring& errorMsg);
    wruntime_error(const wruntime_error& rhs);
    wruntime_error& operator=(const wruntime_error& rhs);
#if defined(WIN32) || defined(_WIN32)
    virtual ~wruntime_error();
#elif defined(__APPLE__)
    virtual ~wruntime_error() throw();
#else
    virtual ~wruntime_error() _GLIBCXX_USE_NOEXCEPT;
#endif
    const std::wstring& errorMsg() const;

  private:
    std::wstring mErrorMsg;
  };
};  // std namespace
