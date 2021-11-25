#ifndef SIMPROP_INTERACTIONS_PSBPHOTODISINTEGRATION_H
#define SIMPROP_INTERACTIONS_PSBPHOTODISINTEGRATION_H

#include <map>
#include <numeric>
#include <string>

#include "simprop/Units.h"
#include "simprop/interactions/AbstractInteraction.h"
#include "simprop/utils/misc.h"

namespace simprop {
namespace interactions {

struct PsbParams {
  double eth_1;
  double eth_2;
  double e0_1;
  double e0_2;
  double xi_1;
  double xi_2;
  double Delta_1;
  double Delta_2;
  double zeta;
};

class PsbPhotoDisintegration : AbstractInteration {
 protected:
  void loadPsbParams(std::string filename);

 public:
  PsbPhotoDisintegration() {
    if (!utils::fileExists(m_tableFilename)) throw std::runtime_error("data file not found");
    loadPsbParams(m_tableFilename);
  }
  virtual ~PsbPhotoDisintegration() = default;

  std::array<double, 3> getCrossSections(PID pid, double photonEnergy) const;

  double getSigma(PID pid, double photonEnergy) const override {
    auto sigmas = getCrossSections(pid, photonEnergy);
    return std::accumulate(sigmas.begin(), sigmas.end(), 0.);
  };

 protected:
  std::string m_tableFilename = "data/xsec_params_pd_PSB.txt";
  std::map<PID, PsbParams> m_psbParams;
  const double m_lowThreshold = 30. * SI::MeV;
  const double m_maxEnergy = 150. * SI::MeV;
  size_t m_nNuclearSpecies = 51;
};

}  // namespace interactions
}  // namespace simprop

#endif