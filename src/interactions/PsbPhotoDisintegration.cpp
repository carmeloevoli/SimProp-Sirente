#include "simprop/interactions/PsbPhotoDisintegration.h"

namespace simprop {
namespace interactions {

const double photodisMinEnergy = 30. * SI::MeV;
// const double e_1 = 30.;     // Photodisintegration minimum energy (MeV)
const double photodisdMaxEnergy = 150. * SI::MeV;
// const double e_max = 150.;  // Photodisintegration maximum energy (MeV)
size_t nNuclearSpecies = 51;

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
    auto csi_1 = p[6];  // adimensional
    auto csi_2 = p[7];  // adimensional
    auto delta_1 = p[8] * SI::MeV;
    auto delta_2 = p[9] * SI::MeV;
    auto zita = p[10];  // adimensional
    auto params = PsbParams{eth_1, eth_2, e0_1, e0_2, csi_1, csi_2, delta_1, delta_2, zita};
    m_psbParams[pid] = params;
  }
  assert(m_psbParams.size() == nNuclearSpecies);
}

double W_function(double eth, double e0, double d0) {
  // here we corrected a typo in eq. (3) of the Stecker & Salamon paper
  const auto e_1 = photodisMinEnergy;
  const auto x1 = (e_1 - e0) / (d0 / M_SQRT2);
  const auto x2 = (e0 - eth) / (d0 / M_SQRT2);
  const auto s_erfx1 = std::erf(x1), s_erfx2 = std::erf(x2);
  return d0 * std::sqrt(M_PI / 8.) * (s_erfx1 + s_erfx2);
}

auto Gaussian = [](double x) { return std::exp(-2. * x * x); };

std::array<double, 3> PsbPhotoDisintegration::getCrossSections(PID pid, double photonEnergy) {
  auto Z = getNucleusChargeNumber(pid);
  auto A = getNucleusMassNumber(pid);
  auto p = m_psbParams[pid];

  double sigma_1 = 0, sigma_2 = 0, sigma_3 = 0;

  const double Sigma_d = 60. * SI::mbarn * SI::MeV * (double)((A - Z) * Z) / (double)(A);

  if ((photonEnergy >= p.eth_1) && (photonEnergy <= photodisMinEnergy)) {
    const double W = W_function(p.eth_1, p.e0_1, p.delta_1);
    sigma_1 = (p.csi_1 * Sigma_d / W) * Gaussian((photonEnergy - p.e0_1) / p.delta_1);
  }

  if ((photonEnergy >= p.eth_2) && (photonEnergy <= photodisMinEnergy)) {
    const double W = W_function(p.eth_2, p.e0_2, p.delta_2);
    sigma_2 = (p.csi_2 * Sigma_d / W) * Gaussian((photonEnergy - p.e0_2) / p.delta_2);
  }

  if ((photonEnergy >= photodisMinEnergy) && (photonEnergy <= photodisdMaxEnergy)) {
    sigma_3 = p.zita * Sigma_d / (photodisdMaxEnergy - photodisMinEnergy);
  }

  return {sigma_1, sigma_2, sigma_3};
}

// std::array<double, 3> PsbPhotoDisintegration::getCrossSections(PID pid, double photonEnergy) {
//   if (photonEnergy >= photodisdMaxEnergy) return {0., 0., 0.};

//   if (Z >= 2 && Z <= 4) {
//     auto params = m_sigmaParams[pid];
//     return getSigmaLowMassNuclei(Z, A, photonEnergy, params);
//   }

//   return {0, 0, 0};
// }

// double Stecker1999PhotoDisintegration::getTotalCrossSections(PID pid, double photonEnergy) {
//   auto s = getCrossSections(pid, photonEnergy);
//   return s[0] + s[1] + s[2];
// }

}  // namespace interactions
}  // namespace simprop