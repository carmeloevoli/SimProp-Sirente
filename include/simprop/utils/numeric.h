#ifndef SIMPROP_UTILS_GSL_H
#define SIMPROP_UTILS_GSL_H

#include <gsl/gsl_integration.h>

#include <cassert>
#include <functional>

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

// Axis
std::vector<double> LinAxis(const double &min, const double &max, const size_t &size);

std::vector<double> LogAxis(const double &min, const double &max, const size_t &size);

double interpolate(double x, const std::vector<double> &X, const std::vector<double> &Y);

double cspline(double x, const std::vector<double> &X, const std::vector<double> &Y);

double interpolateEquidistant(double x, double lo, double hi, const std::vector<double> &Y);

double interpolate2d(double x, double y, const std::vector<double> &X, const std::vector<double> &Y,
                     const std::vector<double> &Z);

template <typename T>
T QAGIntegration(std::function<T(T)> f, T start, T stop, int LIMIT, double rel_error = 1e-4) {
  double a = static_cast<double>(start);
  double b = static_cast<double>(stop);
  double abs_error = 0.0;  // disabled
  int key = GSL_INTEG_GAUSS31;
  double result;
  double error;

  gsl_function F;
  F.function = [](double x, void *vf) -> double {
    auto &func = *static_cast<std::function<double(double)> *>(vf);
    return func(x);
  };
  F.params = &f;

  gsl_integration_workspace *workspace_ptr = gsl_integration_workspace_alloc(LIMIT);
  gsl_integration_qag(&F, a, b, abs_error, rel_error, LIMIT, key, workspace_ptr, &result, &error);
  gsl_integration_workspace_free(workspace_ptr);

  return T(result);
}

template <typename T>
T simpsonIntegration(std::function<T(T)> f, T start, T stop, int N = 100) {
  const T a = start;
  const T b = stop;

  const T h = (b - a) / N;
  const T XI0 = f(a) + f(b);

  T XI1 = 0, XI2 = 0;

  for (int i = 1; i < N; ++i) {
    const T X = a + i * h;
    if (i % 2 == 0)
      XI2 = XI2 + f(X);
    else
      XI1 = XI1 + f(X);
  }

  return h * (XI0 + 2 * XI2 + 4 * XI1) / 3.0;
}

}  // namespace utils
}  // namespace simprop

#endif