#include "simprop/utils/misc.h"

#include <cmath>
#include <fstream>
#include <sstream>
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

size_t countFileLines(const std::string& filename) {
  size_t count = 0;
  std::string line;
  std::ifstream file(filename.c_str());
  while (getline(file, line)) count++;
  return count;
}

bool fileExists(const std::string& filename) {
  std::ifstream f(filename.c_str());
  return f.good();
}

std::vector<std::string> split(std::string s, std::string delimiter) {
  std::vector<std::string> result;
  std::istringstream iss(s);
  for (std::string s; iss >> s;) result.push_back(s);
  return result;
}

}  // namespace utils
}  // namespace simprop