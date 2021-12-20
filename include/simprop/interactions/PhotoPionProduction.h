#ifndef SIMPROP_PHOTOPION_H
#define SIMPROP_PHOTOPION_H

#include <string>

#include "simprop/interactions/AbstractInteraction.h"
#include "simprop/units.h"
#include "simprop/utils/lookupTable.h"
#include "simprop/utils/numeric.h"

namespace simprop {
namespace interactions {

class PhotoPionProduction : AbstractInteration {
 protected:
  const std::string pionSigmaFilename = "data/xsec_ppp.txt";
  utils::LookupTable<45, 1> pionSigma{pionSigmaFilename};

 public:
  PhotoPionProduction() {}
  virtual ~PhotoPionProduction() = default;

  inline double getPhotonEnergyThresholdCoM() const {
    return SI::pionMassC2 + pow2(SI::pionMassC2) / (2 * SI::protonMassC2);
  }

  double getAtS(double s) const;
  double getAtPhotonEnergyCoM(double epsPrime) const;
  double getSigma(PID pid, double photonEnergy) const override;
};

}  // namespace interactions
}  // namespace simprop

#endif