#ifndef INCLUDE_SIMPROP_GIT_REVISION_H
#define INCLUDE_SIMPROP_GIT_REVISION_H

#include <string>

std::string git_sha1();
std::string get_version();
bool git_has_local_changes();

#endif  // INCLUDE_SIMPROP_GIT_REVISION_H