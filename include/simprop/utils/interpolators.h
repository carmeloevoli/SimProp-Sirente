#ifndef SIMPROP_UTILS_INTERPOLATORS_H
#define SIMPROP_UTILS_INTERPOLATORS_H

#include <vector>

namespace simprop {
namespace utils {

double interpolate(double x, const std::vector<double> &X, const std::vector<double> &Y);

double cspline(double x, const std::vector<double> &X, const std::vector<double> &Y);

double interpolateEquidistant(double x, double lo, double hi, const std::vector<double> &Y);

double interpolate2d(double x, double y, const std::vector<double> &X, const std::vector<double> &Y,
                     const std::vector<double> &Z);

}  // namespace utils
}  // namespace simprop

#endif