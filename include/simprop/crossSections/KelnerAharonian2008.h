#ifndef SIMPROP_XSECS_KELNERAHARONIAN2008_H
#define SIMPROP_XSECS_KELNERAHARONIAN2008_H

namespace simprop {
namespace KelnerAharonian2008 {

class SecondarySpectrum {
 public:
  virtual ~SecondarySpectrum() = default;

  double Phi(double eta, double x) {
    const auto z = eta / m_eta_0;
    const auto _B = B(z);
    const auto _s = s(z);
    const auto _delta = delta(z);

    const auto _psi = psi(z);
    const auto _xPrimeMinus = xPrimeMinus(z);
    const auto _xPrimePlus = xPrimePlus(z);

    /* Eq. 32 */
    const auto yPrime = (x - _xPrimeMinus) / (_xPrimePlus - _xPrimeMinus);

    double value = 0.;
    if (x < _xPrimeMinus) {
      /* Eq. 33 */
      value = _B * std::pow(M_LN2, _psi);
    } else if (x < _xPrimePlus) {
      /* Eq. 31 */
      value = _B;
      value *= std::exp(-_s * std::pow(std::log(x / _xPrimeMinus), _delta));
      value *= std::pow(std::log(2. / (1. + pow2(yPrime))), _psi);
    }
    return value;
  }

 protected:
  virtual double B(double z) const = 0;
  virtual double s(double z) const = 0;
  virtual double delta(double z) const = 0;
  virtual double psi(double z) const = 0;
  virtual double xPrimeMinus(double z) const = 0;
  virtual double xPrimePlus(double z) const = 0;

 protected:
  const double m_r = SI::pionMassC2 / SI::protonMassC2;
  const double m_r2 = m_r * m_r;
  /* Eq. 16 */
  const double m_eta_0 = 2. * m_r + m_r2;
};

class AntiNuMuSecondarySpectrum : public SecondarySpectrum {
 public:
  AntiNuMuSecondarySpectrum() {
    std::string filename = "data/xsecs_KA2018_antiNuMu.txt";
    m_s_table.loadTable(filename, 1);
    m_delta_table.loadTable(filename, 2);
    m_B_table.loadTable(filename, 3);
  }

 protected:
  double B(double z) const override {
    double value = 0;
    if (m_B_table.xIsInside(z)) value = m_B_table.get(z);
    return value * SI::cm3 / SI::sec;
  }

  double s(double z) const override {
    double value = 0;
    if (m_s_table.xIsInside(z)) value = m_s_table.get(z);
    return value;
  }

  double delta(double z) const override {
    double value = 0;
    if (m_delta_table.xIsInside(z)) value = m_delta_table.get(z);
    return value;
  }

  double psi(double z) const override {
    /* Eq. 34 */
    return 2.5 + 1.4 * std::log(z);
  }

  double xPrimeMinus(double eta) const override {
    /* Eq. 19 */
    auto xMinus = 1. / 2. / (1. + eta);
    xMinus *= eta + m_r2 - std::sqrt((eta - m_r2 - 2. * m_r) * (eta - m_r2 + 2. * m_r));
    /* Eq. 35 */
    return xMinus / 4.;
  }

  double xPrimePlus(double eta) const override {
    /* Eq. 19 */
    auto xPlus = 1. / 2. / (1. + eta);
    xPlus *= eta + m_r2 + std::sqrt((eta - m_r2 - 2. * m_r) * (eta - m_r2 + 2. * m_r));
    /* Eq. 35 */
    return xPlus;
  }

 protected:
  utils::LookupArray<20> m_B_table;
  utils::LookupArray<20> m_s_table;
  utils::LookupArray<20> m_delta_table;
};

}  // namespace KelnerAharonian2008
}  // namespace simprop

#endif