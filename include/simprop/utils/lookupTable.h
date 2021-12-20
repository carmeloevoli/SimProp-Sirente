#ifndef SIMPROP_UTILS_LOOKUPTABLE_H
#define SIMPROP_UTILS_LOOKUPTABLE_H

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "simprop/utils/io.h"
#include "simprop/utils/numeric.h"

namespace simprop {
namespace utils {

typedef std::vector<double>::const_iterator const_iterator;

template <size_t xSize, size_t ySize>
class LookupTable {
 public:
  explicit LookupTable(std::string filePath) : m_filePath(filePath) {
    if (!utils::fileExists(filePath))
      throw std::runtime_error("file data for lookup table does not exist");
    if (xSize < 2) throw std::runtime_error("x-axis size must be > 1");
    if (ySize < 1) throw std::runtime_error("y-axis size must be > 0");
    m_xAxis.reserve(xSize);
    m_yAxis.reserve(ySize);
    m_table.reserve(xSize * ySize);
    loadTable();
  }

  double get(double x, double y) const {
    return utils::interpolate2d(x, y, m_xAxis, m_yAxis, m_table);
  }
  double get(double x) const { return utils::interpolate(x, m_xAxis, m_table); }
  double spline(double x) const { return utils::cspline(x, m_xAxis, m_table); }

  bool isWithinXRange(double x) const { return x >= m_xAxis.front() && x <= m_xAxis.back(); }
  bool isWithinYRange(double y) const { return y >= m_yAxis.front() && y <= m_yAxis.back(); }

 protected:
  void loadTable() {
    if (ySize == 1)
      loadTable1D();
    else
      loadTable2D();
  }

  void loadTable1D() {
    auto v = utils::loadFileByRow(m_filePath, ",");
    for (size_t i = 0; i < xSize; ++i) {
      auto line = v.at(i);
      if (line.size() != 2) throw std::runtime_error("error in reading table values");
      m_xAxis.emplace_back(line[0]);
      m_table.emplace_back(line[1]);
    }
    assert(m_xAxis.size() == xSize && m_table.size() == xSize);
  }

  void loadTable2D() {
    auto v = utils::loadFileByRow(m_filePath, ",");
    size_t counter = 0;
    for (size_t i = 0; i < xSize; ++i) {
      for (size_t j = 0; j < ySize; ++j) {
        auto line = v.at(counter);
        if (line.size() != 3) throw std::runtime_error("error in reading table values");
        if (j == 0) m_xAxis.emplace_back(line[0]);
        if (i == 0) m_yAxis.emplace_back(line[1]);
        m_table.emplace_back(line[2]);
        counter++;
      }
    }
    assert(m_xAxis.size() == xSize && m_yAxis.size() == ySize);
  }

 protected:
  std::vector<double> m_xAxis;
  std::vector<double> m_yAxis;
  std::vector<double> m_table;
  std::string m_filePath;
};

}  // namespace utils
}  // namespace simprop

#endif