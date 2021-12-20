#include "simprop/utils/numeric.h"

#include <algorithm>
#include <cmath>

#include "simprop/utils/spline.h"

namespace simprop {
namespace utils {

#define index(i, j) ((j) + (i)*Y.size())

double interpolate(double x, const std::vector<double> &X, const std::vector<double> &Y) {
  std::vector<double>::const_iterator it = std::upper_bound(X.begin(), X.end(), x);
  if (it == X.begin()) return 0;
  if (it == X.end()) return 0;

  const size_t i = it - X.begin() - 1;
  return Y[i] + (x - X[i]) * (Y[i + 1] - Y[i]) / (X[i + 1] - X[i]);
}

double interpolateEquidistant(double x, double lo, double hi, const std::vector<double> &Y) {
  if (x <= lo) return 0;
  if (x >= hi) return 0;

  const double dx = (hi - lo) / (Y.size() - 1);
  const double p = (x - lo) / dx;
  const size_t i = std::floor(p);
  return Y[i] + (p - i) * (Y[i + 1] - Y[i]);
}

double cspline(double x, const std::vector<double> &X, const std::vector<double> &Y) {
  if (x < X.front()) return 0;
  if (x > X.back()) return 0;
  tk::spline s(X, Y, tk::spline::cspline);
  return s(x);
}

double interpolate2d(double x, double y, const std::vector<double> &X, const std::vector<double> &Y,
                     const std::vector<double> &Z) {
  std::vector<double>::const_iterator itx = std::upper_bound(X.begin(), X.end(), x);
  std::vector<double>::const_iterator ity = std::upper_bound(Y.begin(), Y.end(), y);

  if (x > X.back() || x < X.front()) return 0;
  if (y > Y.back() || y < Y.front()) return 0;

  if (itx == X.begin() && ity == Y.begin()) return Z.front();
  if (itx == X.end() && ity == Y.end()) return Z.back();

  const size_t i = itx - X.begin() - 1;
  const size_t j = ity - Y.begin() - 1;

  const double Q11 = Z[index(i, j)];
  const double Q12 = Z[index(i, j + 1)];
  const double Q21 = Z[index(i + 1, j)];
  const double Q22 = Z[index(i + 1, j + 1)];

  const double R1 =
      ((X[i + 1] - x) / (X[i + 1] - X[i])) * Q11 + ((x - X[i]) / (X[i + 1] - X[i])) * Q21;
  const double R2 =
      ((X[i + 1] - x) / (X[i + 1] - X[i])) * Q12 + ((x - X[i]) / (X[i + 1] - X[i])) * Q22;

  return ((Y[j + 1] - y) / (Y[j + 1] - Y[j])) * R1 + ((y - Y[j]) / (Y[j + 1] - Y[j])) * R2;
}

}  // namespace utils
}  // namespace simprop
