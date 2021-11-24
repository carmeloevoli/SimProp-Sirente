#ifndef SIMPROP_INTERACTIONS_PSBPHOTODISINTEGRATION_H
#define SIMPROP_INTERACTIONS_PSBPHOTODISINTEGRATION_H

#include <map>
#include <string>

#include "simprop/Units.h"
#include "simprop/interactions/AbstractPhotoDisintegration.h"
#include "simprop/utils/misc.h"

namespace simprop {
namespace interactions {

struct PsbParams {
  double eth_1;
  double eth_2;
  double e0_1;
  double e0_2;
  double csi_1;
  double csi_2;
  double delta_1;
  double delta_2;
  double zita;
};

class PsbPhotoDisintegration : AbstractPhotoDisintegration {
 protected:
  void loadPsbParams(std::string filename);

 public:
  PsbPhotoDisintegration() {
    if (!utils::fileExists(m_tableFilename)) throw std::runtime_error("data file not found");
    loadPsbParams(m_tableFilename);
  }
  virtual ~PsbPhotoDisintegration() = default;

  virtual double getAbsorptionSigma(double photonEnergy) const override { return 0.; };

  std::array<double, 3> getCrossSections(PID pid, double photonEnergy);
  double getTotalCrossSections(PID pid, double photonEnergy);

 protected:
  std::string m_tableFilename = "data/xsec_params_pd_PSB.txt";
  std::map<PID, PsbParams> m_psbParams;
};

}  // namespace interactions
}  // namespace simprop

#endif