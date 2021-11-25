#include "simprop/interactions/PsbPhotoDisintegration.h"

namespace simprop {
namespace interactions {

void PsbPhotoDisintegration::loadPsbParams(std::string filename) {
  auto table = utils::loadFileByRow(filename, ",");
  for (auto &p : table) {
    auto Z = (int)p[0];
    auto A = (int)p[1];
    auto pid = getPidNucleus(Z, A);
    auto eth_1 = p[2] * SI::MeV;
    auto eth_2 = p[3] * SI::MeV;
    auto e0_1 = p[4] * SI::MeV;
    auto e0_2 = p[5] * SI::MeV;
    auto xi_1 = p[6];  // adimensional
    auto xi_2 = p[7];  // adimensional
    auto Delta_1 = p[8] * SI::MeV;
    auto Delta_2 = p[9] * SI::MeV;
    auto zeta = p[10];  // adimensional
    auto params = PsbParams{eth_1, eth_2, e0_1, e0_2, xi_1, xi_2, Delta_1, Delta_2, zeta};
    m_psbParams[pid] = params;
  }
  assert(m_psbParams.size() == m_nNuclearSpecies);
}

double W_function(double e_min, double e_0, double e_1, double Delta) {
  assert(e_min < e_0 && e_0 < e_1);
  const auto x1 = M_SQRT2 * (e_1 - e_0) / Delta;
  const auto x2 = M_SQRT2 * (e_0 - e_min) / Delta;
  const auto Erf_x1 = std::erf(x1), Erf_x2 = std::erf(x2);
  return Delta * std::sqrt(M_PI / 8.) * (Erf_x1 + Erf_x2);
}

double lowEnergyModel(double epsPrime, double xi, double e_0, double e_1, double e_min,
                      double Delta) {
  const double W = W_function(e_min, e_0, e_1, Delta);
  const double t = (epsPrime - e_0) / Delta;
  return xi / W * std::exp(-2. * (t * t));
}

std::array<double, 3> PsbPhotoDisintegration::getCrossSections(
    PID pid, double photonEnergy) const {  // TODO specify frame
  auto Z = getNucleusChargeNumber(pid);
  auto A = getNucleusMassNumber(pid);
  auto p = m_psbParams.at(pid);

  double sSingleNucleon = 0, sDoubleNucleon = 0, sMultipleNucleon = 0;

  if (photonEnergy <= m_maxEnergy) {
    const double Sigma_d = 60. * SI::mbarn * SI::MeV * (double)((A - Z) * Z) / (double)(A);

    if (photonEnergy >= m_lowThreshold) {
      sMultipleNucleon = p.zeta * Sigma_d / (m_maxEnergy - m_lowThreshold);
    } else {
      if (photonEnergy >= p.eth_1) {
        sSingleNucleon = Sigma_d * lowEnergyModel(photonEnergy, p.xi_1, p.e0_1, m_lowThreshold,
                                                  p.eth_1, p.Delta_1);
      }
      if (photonEnergy >= p.eth_2) {
        sDoubleNucleon = Sigma_d * lowEnergyModel(photonEnergy, p.xi_2, p.e0_2, m_lowThreshold,
                                                  p.eth_2, p.Delta_2);
      }
    }
  }
  return {sSingleNucleon, sDoubleNucleon, sMultipleNucleon};
}

}  // namespace interactions
}  // namespace simprop