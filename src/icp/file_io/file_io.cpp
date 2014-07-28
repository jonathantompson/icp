#include <sys\stat.h>
#include <sstream>
#include <string>
#include "icp/file_io/file_io.h"
#include "icp/math/math_types.h"  // for uint

namespace icp {
namespace file_io {
  bool fileExists(const std::string& filename) {
    // TODO: This is a pretty stupid way to check if a file exists.  I think
    // opening a file handler is probably slow.  Rethink this.
    std::ifstream file(filename.c_str(), std::ios::in | std::ios::binary);
    bool ret_val = false;
    if (file.is_open()) {
      ret_val = true;
      file.close();
    }
    return ret_val;
  }

  PathType getPathType(const std::string& path) {
    struct stat s;
    if (stat(path.c_str(), &s) == 0) {
      if (s.st_mode & S_IFDIR) {
        return DIRECTORY_PATH;
      }
      else if(s.st_mode & S_IFREG) {
        return FILE_PATH;
      } else {
        return UNKNOWN_PATH;  // It's something else
      }
    } else {
      return UNKNOWN_PATH;  // Something went wrong
    }
  }

}  // namespace file_io
}  // namespace icp
