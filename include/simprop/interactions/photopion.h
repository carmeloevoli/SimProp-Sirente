#ifndef SIMPROP_PHOTOPION_H
#define SIMPROP_PHOTOPION_H

#include <string>

#include "simprop/utils/lookupTable.h"
#include "simprop/utils/misc.h"

namespace simprop {
namespace interactions {

class PhotoPion {
 protected:
  const std::string pionLowSigmaFilename = "data/sigma_photopion_SOFIA.txt";
  utils::LookupTable<38, 1> pionSigmaLowS{pionLowSigmaFilename};

  const std::string pionHighSigmaFilename = "data/sigma_photopion_highs.txt";
  utils::LookupTable<7, 1> pionSigmaHighS{pionHighSigmaFilename};

 public:
  PhotoPion() {}
  virtual ~PhotoPion() = default;

  double sigma(double E) {
    // E = photon energy in NRF (MeV)
    if (E < SI::pionProductionThreshold) return 0;
    double value = 0;
    const auto s = utils::pow<2>(SI::protonMassC2) + 2 * SI::protonMassC2 * E;
    if (s <= 5. * SI::GeV2) {
      value = pionSigmaLowS.spline(s / SI::GeV2);
    } else {
      const auto logsqrts = std::log10(std::sqrt(s / SI::GeV2));
      const auto logsigma = pionSigmaHighS.spline(std::min(logsqrts, 4.54));
      value = (logsigma < 0.) ? std::pow(10., logsigma) : 0.;
    }
    return std::max(value, 0.) * SI::mbarn;
  }
};

}  // namespace interactions
}  // namespace simprop

#endif