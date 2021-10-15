#ifndef INCLUDE_SIMPROP_LOGGING_H
#define INCLUDE_SIMPROP_LOGGING_H

#include <plog/Appenders/ColorConsoleAppender.h>
#include <plog/Appenders/ConsoleAppender.h>
#include <plog/Appenders/RollingFileAppender.h>
#include <plog/Formatters/CsvFormatter.h>
#include <plog/Formatters/TxtFormatter.h>
#include <plog/Init.h>
#include <plog/Log.h>

#include "simprop/git_revision.h"

namespace SimProp {
namespace log {

void startup_information();

}  // namespace log
}  // namespace SimProp

#endif  // INCLUDE_SIMPROP_LOGGING_H