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
namespace log {

void startup_information();

}  // namespace log
}  // namespace simprop

#endif  // SIMPROP_UTILS_LOGGING_H