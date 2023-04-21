#ifndef SIMPROP_XSECS_KELNERAHARONIAN2008_H
#define SIMPROP_XSECS_KELNERAHARONIAN2008_H

#include <string>
#include <vector>

#include "simprop/core/units.h"

namespace simprop {
namespace KelnerAharonian2008 {

class SecondarySpectrum {
 public:
  SecondarySpectrum() {}
  virtual ~SecondarySpectrum() = default;
  double Phi(double eta, double x) const;

 protected:
  double xMinus(double eta) const;
  double xPlus(double eta) const;
  double B(double rho) const;
  double s(double rho) const;
  double delta(double rho) const;
  void loadTables(std::string filename);

 protected:
  virtual double psi(double rho) const = 0;
  virtual double xPrimeMinus(double eta) const = 0;
  virtual double xPrimePlus(double eta) const = 0;

 protected:
  const double m_r = SI::pionMassC2 / SI::protonMassC2;
  const double m_r2 = m_r * m_r;
  const double m_eta_0 = 2. * m_r + m_r2;

 protected:
  std::vector<double> m_rho_table;
  std::vector<double> m_s_table;
  std::vector<double> m_delta_table;
  std::vector<double> m_B_table;
};

class AntiNuMuSpectrum final : public SecondarySpectrum {
 public:
  AntiNuMuSpectrum();

 protected:
  double psi(double z) const override;
  double xPrimeMinus(double eta) const override;
  double xPrimePlus(double eta) const override;
};

class NuMuSpectrum final : public SecondarySpectrum {
 public:
  NuMuSpectrum();

 protected:
  double psi(double z) const override;
  double xPrimeMinus(double eta) const override;
  double xPrimePlus(double eta) const override;
};

class AntiNuElectronSpectrum final : public SecondarySpectrum {
 public:
  AntiNuElectronSpectrum();

 protected:
  double psi(double z) const override;
  double xPrimeMinus(double eta) const override;
  double xPrimePlus(double eta) const override;
};

class NuElectronSpectrum final : public SecondarySpectrum {
 public:
  NuElectronSpectrum();

 protected:
  double psi(double z) const override;
  double xPrimeMinus(double eta) const override;
  double xPrimePlus(double eta) const override;
};

struct NeutrinoProductionSpectrum {
  NuMuSpectrum numu;
  AntiNuMuSpectrum antiNumu;
  NuElectronSpectrum nue;
  AntiNuElectronSpectrum antiNue;

  double Phi(double eta, double x) {
    return numu.Phi(eta, x) + antiNumu.Phi(eta, x) + nue.Phi(eta, x) + antiNue.Phi(eta, x);
  }
};

}  // namespace KelnerAharonian2008
}  // namespace simprop

#endif