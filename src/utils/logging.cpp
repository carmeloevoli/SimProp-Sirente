#include "simprop/utils/logging.h"

namespace simprop {
namespace log {

void startup_information() {
  static plog::RollingFileAppender<plog::CsvFormatter> fileAppender("simproplog.csv", 8000, 3);
  static plog::ConsoleAppender<plog::TxtFormatter> consoleAppender;
#ifdef DEBUG
  plog::init(plog::debug, &fileAppender).addAppender(&consoleAppender);
#else
  plog::init(plog::info, &fileAppender).addAppender(&consoleAppender);
#endif
  LOGI << "Welcome to SimProp version " << get_version();
  LOGI << "was built on " << __DATE__ << " at " << __TIME__;
  LOGI << "git version is " << git_sha1();
  LOGW << "has local changes " << git_has_local_changes();
}

}  // namespace log
}  // namespace simprop