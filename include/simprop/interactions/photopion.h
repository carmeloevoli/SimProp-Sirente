#ifndef SIMPROP_PHOTOPION_H
#define SIMPROP_PHOTOPION_H

#include <string>

#include "simprop/utils/lookupTable.h"
#include "simprop/utils/misc.h"

namespace simprop {
namespace interactions {

static constexpr double pionProductionThreshold =
    SI::pionMassC2 + pow2(SI::pionMassC2) / (2 * SI::protonMassC2);

class PhotoPion {
 protected:
  const std::string pionLowSigmaFilename = "data/sigma_photopion_SOFIA.txt";
  utils::LookupTable<38, 1> pionSigmaLowS{pionLowSigmaFilename};

  const std::string pionHighSigmaFilename = "data/sigma_photopion_highs.txt";
  utils::LookupTable<7, 1> pionSigmaHighS{pionHighSigmaFilename};

 public:
  PhotoPion() {}
  virtual ~PhotoPion() = default;

  double sigma_of_s(double s) const {
    double value = 0;
    if (s <= 5. * SI::GeV2) {
      value = pionSigmaLowS.get(s / SI::GeV2);
    } else {
      const auto logsqrts = std::log10(std::sqrt(s / SI::GeV2));
      const auto logsigma = pionSigmaHighS.spline(std::min(logsqrts, 4.54));
      value = (logsigma < 0.) ? std::pow(10., logsigma) : 0.;
    }
    return std::max(value, 0.) * SI::mbarn;
  }

  double sigma(double E) const {
    // E = photon energy in NRF (MeV)
    if (E < pionProductionThreshold) return 0;
    const auto s = utils::pow<2>(SI::protonMassC2) + 2 * SI::protonMassC2 * E;
    return sigma_of_s(s);
  }
};

}  // namespace interactions
}  // namespace simprop

#endif