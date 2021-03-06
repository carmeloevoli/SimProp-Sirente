#include "simprop/utils/logging.h"

namespace simprop {
namespace utils {

void startup_information() {
  static plog::ColorConsoleAppender<plog::TxtFormatter> consoleAppender;
#ifdef DEBUG
  plog::init(plog::debug, &consoleAppender);
#else
  plog::init(plog::info, &consoleAppender);
#endif

  LOGI << "Welcome to SimProp version " << get_version();
  LOGI << "was built on " << __DATE__ << " at " << __TIME__;
  LOGI << "git version is " << git_sha1();
  LOGW << "has local changes " << git_has_local_changes();
}

}  // namespace utils
}  // namespace simprop