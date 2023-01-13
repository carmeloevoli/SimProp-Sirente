#ifndef SIMPROP_UTILS_LOOKUPTABLE_H
#define SIMPROP_UTILS_LOOKUPTABLE_H

#include <vector>

#include "simprop/core/units.h"
#include "simprop/utils/io.h"
#include "simprop/utils/numeric.h"
#include "simprop/utils/timer.h"

namespace simprop {
namespace utils {

// typedef std::vector<double>::const_iterator const_iterator;

template <size_t xSize>
class LookupArray {
 public:
  LookupArray() {
    if (xSize < 2) throw std::runtime_error("x-axis size must be > 1");
    m_xAxis.reserve(xSize);
    m_array.reserve(xSize);
  }

  inline double get(double x) const { return utils::interpolate(x, m_xAxis, m_array); }
  inline double spline(double x) const { return utils::cspline(x, m_xAxis, m_array); }
  inline bool xIsInside(double x) const { return x >= m_xAxis.front() && x <= m_xAxis.back(); }

 public:
  void loadTable(const std::string& filePath, size_t iCol = 1) {
    if (!utils::fileExists(filePath))
      throw std::runtime_error("file data for lookup array does not exist");
    auto v = utils::loadFileByRow(filePath, ",");
    for (size_t i = 0; i < xSize; ++i) {
      auto line = v.at(i);
      m_xAxis.emplace_back(line[0]);
      m_array.emplace_back(line[iCol]);
    }
    assert(m_xAxis.size() == xSize && m_array.size() == xSize);
  }

  void cacheTable(const std::function<double(double)>& func,
                  const std::pair<double, double>& range) {
    const double dx = (range.second - range.first) / (double)(xSize - 1);
    utils::Timer timer("time for caching");
    for (size_t i = 0; i < xSize; ++i) {
      auto x = (double)i * dx + range.first;
      auto f_x = func(x);
      m_xAxis.emplace_back(x);
      m_array.emplace_back(f_x);
    }
    assert(m_xAxis.size() == xSize && m_array.size() == xSize);
  }

 protected:
  std::vector<double> m_xAxis;
  std::vector<double> m_array;
};

// template <size_t xSize, size_t ySize>
// class LookupTable {
//  public:
//   explicit LookupTable(std::string filePath) : m_filePath(filePath) {
//     if (!utils::fileExists(filePath))
//       throw std::runtime_error("file data for lookup table does not exist");
//     if (xSize < 2) throw std::runtime_error("x-axis size must be > 1");
//     if (ySize < 2) throw std::runtime_error("y-axis size must be > 1");
//     m_xAxis.reserve(xSize);
//     m_yAxis.reserve(ySize);
//     m_table.reserve(xSize * ySize);
//     loadTable();
//   }

//   double get(double x, double y) const {
//     return utils::interpolate2d(x, y, m_xAxis, m_yAxis, m_table);
//   }

//   bool xIsInside(double x) const { return x >= m_xAxis.front() && x <= m_xAxis.back(); }
//   bool yIsInside(double y) const { return y >= m_yAxis.front() && y <= m_yAxis.back(); }

//  protected:
//   void loadTable() {
//     auto v = utils::loadFileByRow(m_filePath, ",");
//     size_t counter = 0;
//     for (size_t i = 0; i < xSize; ++i) {
//       for (size_t j = 0; j < ySize; ++j) {
//         auto line = v.at(counter);
//         if (line.size() != 3) throw std::runtime_error("error in reading table values");
//         if (j == 0) m_xAxis.emplace_back(line[0]);
//         if (i == 0) m_yAxis.emplace_back(line[1]);
//         m_table.emplace_back(line[2]);
//         counter++;
//       }
//     }
//     assert(m_xAxis.size() == xSize && m_yAxis.size() == ySize);
//   }

//  protected:
//   std::vector<double> m_xAxis;
//   std::vector<double> m_yAxis;
//   std::vector<double> m_table;
//   std::string m_filePath;
// };

}  // namespace utils
}  // namespace simprop

#endif