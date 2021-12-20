#ifndef SIMPROP_UTILS_LOGGING_H
#define SIMPROP_UTILS_LOGGING_H

#include <plog/Appenders/ColorConsoleAppender.h>
#include <plog/Appenders/ConsoleAppender.h>
#include <plog/Appenders/RollingFileAppender.h>
#include <plog/Formatters/CsvFormatter.h>
#include <plog/Formatters/TxtFormatter.h>
#include <plog/Init.h>
#include <plog/Log.h>

#include "simprop/utils/git_revision.h"

namespace simprop {
namespace utils {

void startup_information();

}  // namespace utils
}  // namespace simprop

#endif  // SIMPROP_UTILS_LOGGING_H