#include "simprop/utils.h"

#include <stdexcept>

namespace simprop {
namespace utils {

std::string removeExtensionIniFilename(std::string inputFilename) {
  auto ext = inputFilename.substr(inputFilename.rfind('.') + 1);
  if (ext != "ini") throw std::runtime_error("Wrong input filename! It must be : <filename>.ini");
  size_t lastindex = inputFilename.find_last_of(".");
  return inputFilename.substr(0, lastindex);
}

}  // namespace utils
}  // namespace simprop