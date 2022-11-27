#include "simprop/utils/io.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "simprop/utils/logging.h"

namespace simprop {
namespace utils {

std::string removeExtensionIniFilename(std::string inputFilename) {
  auto ext = inputFilename.substr(inputFilename.rfind('.') + 1);
  if (ext != "ini")
    throw std::invalid_argument("Wrong input filename! It must be : <filename>.ini");
  size_t lastindex = inputFilename.find_last_of(".");
  return inputFilename.substr(0, lastindex);
}

size_t countFileLines(const std::string& filename) {
  size_t counter = 0;
  std::string line;
  std::ifstream file(filename.c_str());
  while (getline(file, line)) {
    if (line.at(0) != '#') counter++;
  }
  return counter;
}

bool fileExists(const std::string& filename) {
  std::ifstream f(filename.c_str());
  return f.good();
}

std::vector<std::string> split(std::string s, std::string delimiter) {
  std::vector<std::string> v;
  size_t pos = 0;
  std::string token;
  while ((pos = s.find(delimiter)) != std::string::npos) {
    token = s.substr(0, pos);
    v.push_back(token);
    s.erase(0, pos + delimiter.length());
  }
  v.push_back(s);
  return v;
}

std::vector<double> loadRow(std::string filePath, size_t iRow, std::string delimiter) {
  if (iRow > countFileLines(filePath)) throw std::runtime_error("row index outside file size");
  std::vector<double> v;
  size_t count = 0;
  std::string line;
  std::ifstream file(filePath.c_str());
  while (getline(file, line)) {
    if (line.at(0) != '#') {
      if (iRow == count) {
        auto s = split(line, delimiter);
        v.resize(s.size());
        std::transform(s.begin(), s.end(), v.begin(),
                       [](const std::string& value) { return std::stod(value); });
      }
      count++;
    }
  }
  return v;
}

std::vector<std::vector<double> > loadFileByRow(std::string filePath, std::string delimiter) {
  std::vector<std::vector<double> > rows;
  size_t count = 0;
  std::string line;
  std::ifstream file(filePath.c_str());
  while (getline(file, line)) {
    if (line.at(0) != '#') {
      std::vector<double> v;
      auto s = split(line, delimiter);
      v.resize(s.size());
      std::transform(s.begin(), s.end(), v.begin(),
                     [](const std::string& value) { return std::stod(value); });
      rows.push_back(v);
    }
    count++;
  }
  return rows;
}

OutputFile::OutputFile(const std::string& name) : filename(name), out("output/" + name) {}
OutputFile::~OutputFile() { LOGD << "created output file " << filename; }

}  // namespace utils
}  // namespace simprop