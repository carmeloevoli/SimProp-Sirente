#ifndef SIMPROP_UTILS_LOOKUPTABLE_H
#define SIMPROP_UTILS_LOOKUPTABLE_H

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "simprop/utils/interpolators.h"
#include "simprop/utils/misc.h"

namespace simprop {
namespace utils {

template <size_t xSize, size_t ySize>
class LookupTable {
 public:
  explicit LookupTable(std::string filePath) : m_filePath(filePath) {
    if (!utils::fileExists(filePath))
      throw std::runtime_error("file data for lookup table does not exist");
    if (xSize < 2) throw std::runtime_error("x-axis size must be > 1");
    loadXAxis();
    if (ySize > 1) loadYAxis();
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
  void loadXAxis() {
    auto v = utils::loadRow(m_filePath, 0, ",");
    if (v.size() != xSize) throw std::runtime_error("error in reading x-axis");
    if (!std::is_sorted(v.begin(), v.end())) throw std::runtime_error("x-axis in not sorted");
    m_xAxis.reserve(xSize);
    std::copy(v.begin(), v.end(), std::back_inserter(m_xAxis));
  }
  void loadYAxis() {
    auto v = utils::loadRow(m_filePath, 1, ",");
    if (v.size() != ySize) throw std::runtime_error("error in reading y-axis");
    m_yAxis.reserve(ySize);
    std::copy(v.begin(), v.end(), std::back_inserter(m_yAxis));
  }
  void loadTable() {
    m_table.reserve(xSize * ySize);
    for (size_t i = 0; i < ySize; ++i) {
      auto v = utils::loadRow(m_filePath, i + 1 + (ySize > 1), ",");
      if (v.size() != xSize) throw std::runtime_error("error in reading table values");
      std::copy(v.begin(), v.end(), std::back_inserter(m_table));
    }
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