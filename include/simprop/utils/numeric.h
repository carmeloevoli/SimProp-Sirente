#ifndef SIMPROP_UTILS_GSL_H
#define SIMPROP_UTILS_GSL_H

#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_roots.h>

#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
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
T deriv(std::function<T(T)> f, T x, double rel_error = 1e-4) {
  double result;
  double abs_error = 0.0;  // disabled
  gsl_function F;

  F.function = [](double x, void *vf) -> double {
    auto &func = *static_cast<std::function<double(double)> *>(vf);
    return func(x);
  };
  F.params = &f;

  gsl_deriv_central(&F, static_cast<double>(x), rel_error, &result, &abs_error);

  return T(result);
}

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

template <typename T>
T rk4fixed(std::function<T(T, T)> dydx, T yStart, T xStart, T xEnd, T h) {
  size_t steps = (size_t)((xEnd - xStart) / h);
  T x = xStart;
  T y = yStart;
  while (steps--) {
    auto k1 = h * dydx(x, y);
    auto k2 = h * dydx(x + 0.5 * h, y + 0.5 * k1);
    auto k3 = h * dydx(x + 0.5 * h, y + 0.5 * k2);
    auto k4 = h * dydx(x + h, y + k3);
    y += (k1 + 2. * k2 + 2. * k3 + k4) / 6.;
    x += h;
  }
  return y;
}

template <typename T>
T odeiv(std::function<T(T, T)> dydx, T yStart, T xStart, T xEnd, T rel_error = 1e-4) {
  const size_t NEQS = 1;
  double abs_error = 0.0;  // disabled

  const gsl_odeiv2_step_type *stepType = gsl_odeiv2_step_rkf45;
  gsl_odeiv2_step *s = gsl_odeiv2_step_alloc(stepType, NEQS);
  gsl_odeiv2_control *c = gsl_odeiv2_control_y_new(abs_error, static_cast<double>(rel_error));
  gsl_odeiv2_evolve *e = gsl_odeiv2_evolve_alloc(NEQS);

  double y[1] = {yStart};

  auto func = [](double t, const double *y, double *f, void *params) -> int {
    auto dydx = *static_cast<std::function<double(double, double)> *>(params);
    f[0] = dydx(t, y[0]);
    return GSL_SUCCESS;
  };

  gsl_odeiv2_system sys = {func, nullptr, NEQS, &dydx};

  double t = xStart, t1 = xEnd;
  double h = 1e-5 * (xEnd - xStart);

  while (t < t1) {
    int status = gsl_odeiv2_evolve_apply(e, c, s, &sys, &t, t1, &h, y);
    if (status != GSL_SUCCESS) break;
  }

  gsl_odeiv2_evolve_free(e);
  gsl_odeiv2_control_free(c);
  gsl_odeiv2_step_free(s);
  return y[0];
}

}  // namespace utils
}  // namespace simprop

#endif