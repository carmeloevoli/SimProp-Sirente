#include "simprop/utils/misc.h"

#include <cmath>
#include <stdexcept>

namespace simprop {
namespace utils {

std::string removeExtensionIniFilename(std::string inputFilename) {
  auto ext = inputFilename.substr(inputFilename.rfind('.') + 1);
  if (ext != "ini")
    throw std::invalid_argument("Wrong input filename! It must be : <filename>.ini");
  size_t lastindex = inputFilename.find_last_of(".");
  return inputFilename.substr(0, lastindex);
}

std::vector<double> LinAxis(const double& min, const double& max, const size_t& size) {
  if (!(min < max)) throw std::invalid_argument("min must be smaller than max");
  if (!(size > 1)) throw std::invalid_argument("size must be larger than 1");

  const double dx = (max - min) / (double)(size - 1);
  std::vector<double> v(size);
  for (size_t i = 0; i < size; ++i) {
    const auto value = min + dx * i;
    v[i] = value;
  }
  return v;
}

std::vector<double> LogAxis(const double& min, const double& max, const size_t& size) {
  if (!(min < max)) throw std::invalid_argument("min must be smaller than max");
  if (!(size > 1)) throw std::invalid_argument("size must be larger than 1");

  const double delta_log = std::exp(std::log(max / min) / (size - 1));
  std::vector<double> v(size);
  for (size_t i = 0; i < size; ++i) {
    const auto value = std::exp(std::log(min) + (double)i * std::log(delta_log));
    v[i] = value;
  }
  return v;
}

}  // namespace utils
}  // namespace simprop