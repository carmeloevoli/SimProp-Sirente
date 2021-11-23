#ifndef SIMPROP_PHOTOPION_H
#define SIMPROP_PHOTOPION_H

#include <string>

#include "simprop/utils/lookupTable.h"
#include "simprop/utils/misc.h"

namespace simprop {
namespace interactions {

class PhotoPionProduction {
 protected:
  const std::string pionSigmaFilename = "data/xsec_ppp.txt";
  utils::LookupTable<45, 1> pionSigma{pionSigmaFilename};

 public:
  PhotoPionProduction() {}
  virtual ~PhotoPionProduction() = default;

  double getPhotonEnergyThresholdCoM() const {
    return SI::pionMassC2 + pow2(SI::pionMassC2) / (2 * SI::protonMassC2);
  }

  double getAtS(double s) const {
    auto logs = std::log10(s / SI::GeV2);
    auto value = pionSigma.get(logs);
    return std::max(value, 0.) * SI::mbarn;
  }

  double getAtPhotonEnergyCoM(double epsPrime) const {
    if (epsPrime < getPhotonEnergyThresholdCoM()) return 0;
    const auto s = utils::pow<2>(SI::protonMassC2) + 2 * SI::protonMassC2 * epsPrime;
    return getAtS(s);
  }
};

}  // namespace interactions
}  // namespace simprop

#endif