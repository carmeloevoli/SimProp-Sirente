#include "simprop/photonFields/Dominguez2011PhotonField.h"

#include "simprop/Units.h"
#include "simprop/common.h"

namespace simprop {
namespace photonfield {

Dominguez2011PhotonField::Dominguez2011PhotonField() {}

double Dominguez2011PhotonField::getPhotonDensity(double ePhoton, double z) const {
  const auto w = energyToWavelenght(ePhoton);
  const auto logw = std::log10(w / SI::micron);
  if (m_field.isWithinXRange(logw) && m_field.isWithinYRange(z)) {
    const auto value = m_field.get(logw, z);
    const auto power = std::pow(10., value) * SI::nW / SI::m2 / SI::sr;
    return power / ePhoton / ePhoton / SI::cOver4pi;
  } else {
    return 0;
  }
};

}  // namespace photonfield
}  // namespace simprop