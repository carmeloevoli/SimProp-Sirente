#ifndef SIMPROP_UTILS_H
#define SIMPROP_UTILS_H

#include <fstream>
#include <string>
#include <vector>

#include "simprop/utils/logging.h"

namespace simprop {
namespace utils {

// pow implementation as template for integer exponents
template <unsigned int exponent>
inline double pow(double base) {
  return pow<(exponent >> 1)>(base * base) * (((exponent & 1) > 0) ? base : 1);
}

template <>
inline double pow<0>(double base) {
  return 1;
}

// String Utilities
std::string removeExtensionIniFilename(std::string inputFilename);

// Axis
std::vector<double> LinAxis(const double& min, const double& max, const size_t& size);
std::vector<double> LogAxis(const double& min, const double& max, const size_t& size);

// Files
size_t countFileLines(const std::string& filename);
bool fileExists(const std::string& filename);
std::vector<std::string> split(std::string s, std::string delimiter = " ");
std::vector<double> loadRow(std::string filePath, size_t iRow, std::string delimiter = " ");

// Output file
class OutputFile {
 protected:
  std::string filename;
  std::ofstream out;

 public:
  explicit OutputFile(std::string filename) {
    this->filename = filename;
    out.open("output/" + filename);
  }
  ~OutputFile() {
    out.close();
    LOGD << "created output file " << filename;
  }
  std::ofstream& operator()() { return out; }
};

}  // namespace utils
}  // namespace simprop

#endif  // SIMPROP_UTILS_H