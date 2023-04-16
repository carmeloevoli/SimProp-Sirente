#include "simprop/crossSections/KelnerAharonian2008.h"

#include <cmath>

#include "simprop/core/units.h"
#include "simprop/utils/io.h"
#include "simprop/utils/numeric.h"

namespace simprop {
namespace KelnerAharonian2008 {

void SecondarySpectrum::loadTables(std::string filename) {
  if (utils::fileExists(filename)) {
    auto nRows = utils::countFileLines(filename);
    auto v = utils::loadFileByRow(filename, ",");
    size_t counter = 0;
    for (size_t i = 0; i < nRows; ++i) {
      auto line = v.at(counter);
      if (line.size() != 4) throw std::runtime_error("error in reading table values");
      m_rho_table.push_back(line[0]);
      m_s_table.push_back(line[1]);
      m_delta_table.push_back(line[2]);
      m_B_table.push_back(line[3] * SI::cm3 / SI::sec);
      counter++;
    }
    assert(m_rho_table.size() == nRows && m_B_table.size() == nRows);
  } else {
    throw std::runtime_error("error reading from file : " + filename);
  }
}

double SecondarySpectrum::B(double rho) const {
  double value = 0;
  if (rho > m_rho_table.front() && rho < m_rho_table.back()) {
    return utils::interpolate(rho, m_rho_table, m_B_table);
  }
  return value;
}

double SecondarySpectrum::s(double rho) const {
  double value = 0;
  if (rho > m_rho_table.front() && rho < m_rho_table.back()) {
    return utils::interpolate(rho, m_rho_table, m_s_table);
  }
  return value;
}

double SecondarySpectrum::delta(double rho) const {
  double value = 0;
  if (rho > m_rho_table.front() && rho < m_rho_table.back()) {
    return utils::interpolate(rho, m_rho_table, m_delta_table);
  }
  return value;
}

double SecondarySpectrum::Phi(double eta, double x) const {
  const auto rho = eta / m_eta_0;
  const auto _B = B(rho);
  const auto _s = s(rho);
  const auto _delta = delta(rho);
  const auto _psi = psi(rho);

  const auto _xPrimeMinus = xPrimeMinus(eta);
  const auto _xPrimePlus = xPrimePlus(eta);

  const auto yPrime = (x - _xPrimeMinus) / (_xPrimePlus - _xPrimeMinus);

  double value = 0.;
  if (x < _xPrimeMinus) {
    value = _B * std::pow(M_LN2, _psi);
  } else if (x < _xPrimePlus) {
    value = _B;
    value *= std::exp(-_s * std::pow(std::log(x / _xPrimeMinus), _delta));
    value *= std::pow(std::log(2. / (1. + pow2(yPrime))), _psi);
  }
  return value;
}

double SecondarySpectrum::xMinus(double eta) const {
  auto value = 1. / 2. / (1. + eta);
  value *= eta + m_r2 - std::sqrt((eta - m_r2 - 2. * m_r) * (eta - m_r2 + 2. * m_r));
  return value;
}

double SecondarySpectrum::xPlus(double eta) const {
  auto value = 1. / 2. / (1. + eta);
  value *= eta + m_r2 + std::sqrt((eta - m_r2 - 2. * m_r) * (eta - m_r2 + 2. * m_r));
  return value;
}

AntiNuMuSpectrum::AntiNuMuSpectrum() { loadTables("data/xsecs_KA2018_antiNuMu.txt"); }

double AntiNuMuSpectrum::psi(double rho) const { return 2.5 + 1.4 * std::log(rho); }

double AntiNuMuSpectrum::xPrimeMinus(double eta) const { return xMinus(eta) / 4.; }

double AntiNuMuSpectrum::xPrimePlus(double eta) const { return xPlus(eta); }

NuMuSpectrum::NuMuSpectrum() { loadTables("data/xsecs_KA2018_nuMu.txt"); }

double NuMuSpectrum::psi(double rho) const { return 2.5 + 1.4 * std::log(rho); }

double NuMuSpectrum::xPrimeMinus(double eta) const { return 0.427 * xMinus(eta); }

double NuMuSpectrum::xPrimePlus(double eta) const {
  const double rho = eta / m_eta_0;
  if (rho < 2.14) {
    return 0.427 * xPlus(eta);
  } else if (rho < 10.) {
    return (0.427 + 0.0729 * (rho - 2.14)) * xPlus(eta);
  } else {
    return xPlus(eta);
  }
}

NuElectronSpectrum::NuElectronSpectrum() { loadTables("data/xsecs_KA2018_nuElectron.txt"); }

double NuElectronSpectrum::psi(double rho) const { return 2.5 + 1.4 * std::log(rho); }

double NuElectronSpectrum::xPrimeMinus(double eta) const { return xMinus(eta) / 4.; }

double NuElectronSpectrum::xPrimePlus(double eta) const { return xPlus(eta); }

AntiNuElectronSpectrum::AntiNuElectronSpectrum() {
  loadTables("data/xsecs_KA2018_antiNuElectron.txt");
}

double AntiNuElectronSpectrum::psi(double rho) const {
  return (rho > 4.) ? 6. * (1. - std::exp(1.5 * (4. - rho))) : 0.;
}

double AntiNuElectronSpectrum::xPrimeMinus(double eta) const {
  double value = 1. / 2. / (1. + eta);
  value *= eta - 2. * m_r - std::sqrt(eta * (eta - 4. * m_r * (1. + m_r)));
  return value;
}

double AntiNuElectronSpectrum::xPrimePlus(double eta) const {
  double value = 1. / 2. / (1. + eta);
  value *= eta - 2. * m_r + std::sqrt(eta * (eta - 4. * m_r * (1. + m_r)));
  return value;
}

double NeutrinoSpectrum::get(double eta, double x) const {
  return numu.Phi(eta, x) + antiNumu.Phi(eta, x) + nue.Phi(eta, x) + antiNue.Phi(eta, x);
}

}  // namespace KelnerAharonian2008
}  // namespace simprop