#ifndef SIMPROP_UTILS_GSL_H
#define SIMPROP_UTILS_GSL_H

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include <cassert>
#include <functional>
#include <vector>

namespace simprop {
namespace utils {

// Axis
std::vector<double> LinAxis(const double &min, const double &max, const size_t &size);

std::vector<double> LogAxis(const double &min, const double &max, const size_t &size);

inline bool isInside(double x, const std::vector<double> &X) {
  return (x >= X.front() && x <= X.back());
};

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

template <typename T>
T rootFinder(std::function<T(T)> f, T xLower, T xUpper, int maxIter, double relError = 1e-4) {
  int status;
  int iter = 0;
  const gsl_root_fsolver_type *solverType;
  gsl_root_fsolver *solver;

  T r = 0;

  gsl_function F;
  F.function = [](double x, void *vf) -> double {
    auto &func = *static_cast<std::function<double(double)> *>(vf);
    return func(x);
  };
  F.params = &f;

  solverType = gsl_root_fsolver_brent;
  solver = gsl_root_fsolver_alloc(solverType);
  gsl_root_fsolver_set(solver, &F, xLower, xUpper);

  do {
    iter++;
    status = gsl_root_fsolver_iterate(solver);
    r = (T)gsl_root_fsolver_root(solver);
    xLower = gsl_root_fsolver_x_lower(solver);
    xUpper = gsl_root_fsolver_x_upper(solver);
    status = gsl_root_test_interval(xLower, xUpper, 0, relError);
  } while (status == GSL_CONTINUE && iter < maxIter);

  gsl_root_fsolver_free(solver);

  return r;
}

}  // namespace utils
}  // namespace simprop

#endif