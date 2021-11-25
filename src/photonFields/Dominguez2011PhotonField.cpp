#include "simprop/photonFields/Dominguez2011PhotonField.h"

#include "simprop/Units.h"
#include "simprop/common.h"
#include "simprop/utils/gsl.h"

namespace simprop {
namespace photonfield {

Dominguez2011PhotonField::Dominguez2011PhotonField()
    : AbstractPhotonField("Dominguez2011", 1.240e-3 * SI::eV, 1.227 * SI::eV) {}

double Dominguez2011PhotonField::getPhotonDensity(double ePhoton, double z) const {
  const auto loge = std::log10(ePhoton / SI::eV);
  if (m_field.isWithinXRange(z) && m_field.isWithinYRange(loge)) {
    const auto value = m_field.get(z, loge);
    const auto density = std::pow(10., value) / SI::eV / SI::m3;
    return density;
  } else {
    return 0;
  }
};

double Dominguez2011PhotonField::I_gamma(double ePhoton, double z) const {
  const auto loge = std::log10(ePhoton / SI::eV);
  if (m_field.isWithinXRange(z) && m_field.isWithinYRange(loge)) {
    const auto value = m_Igamma.get(z, loge);
    const auto I = std::pow(10., value) / utils::pow<2>(SI::eV) / utils::pow<2>(SI::m3);
    return I;
  } else {
    return 0;
  }
}

}  // namespace photonfield
}  // namespace simprop