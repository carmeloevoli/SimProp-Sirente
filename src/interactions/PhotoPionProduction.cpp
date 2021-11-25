#include "simprop/interactions/PhotoPionProduction.h"

namespace simprop {
namespace interactions {

double PhotoPionProduction::getAtS(double s) const {
  auto logs = std::log10(s / SI::GeV2);
  auto value = pionSigma.get(logs);
  return std::max(value, 0.) * SI::mbarn;
}

double PhotoPionProduction::getAtPhotonEnergyCoM(double epsPrime) const {
  if (epsPrime < getPhotonEnergyThresholdCoM()) return 0;
  const auto s = utils::pow<2>(SI::protonMassC2) + 2 * SI::protonMassC2 * epsPrime;
  return getAtS(s);
}

double PhotoPionProduction::getSigma(PID pid, double photonEnergy) const {
  return getAtPhotonEnergyCoM(photonEnergy);
}

}  // namespace interactions
}  // namespace simprop