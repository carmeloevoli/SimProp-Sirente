#include "simprop/photonFields/CmbPhotonField.h"

#include <cmath>

#include "simprop/Units.h"
#include "simprop/utils/misc.h"

namespace simprop {
namespace photonfield {

BlackbodyPhotonField::BlackbodyPhotonField(std::string fieldName, double blackbodyTemperature) {
  m_blackbodyTemperature = blackbodyTemperature;
  m_fieldName = fieldName;
}

double BlackbodyPhotonField::getPhotonDensity(double ePhoton, double z) const {
  auto density = 8 * M_PI * utils::pow<3>(ePhoton / (SI::hPlanck * SI::cLight)) /
                 std::expm1(ePhoton / SI::kBoltzmann / m_blackbodyTemperature / (1. + z));
  return density / ePhoton;
}

}  // namespace photonfield
}  // namespace simprop