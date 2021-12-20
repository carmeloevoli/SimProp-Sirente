#ifndef SIMPROP_PHOTONFIELDS_DOMINGUEZ2011_H
#define SIMPROP_PHOTONFIELDS_DOMINGUEZ2011_H

#include <string>
#include <vector>

#include "simprop/photonFields/AbstractPhotonField.h"
#include "simprop/units.h"
#include "simprop/utils/lookupTable.h"

namespace simprop {
namespace photonfield {

class Dominguez2011PhotonField : public AbstractPhotonField {
 protected:
  utils::LookupTable<18, 50> m_field{"data/EBL_Dominguez2011.txt"};
  utils::LookupTable<18, 50> m_Igamma{"data/EBL_Igamma_Dominguez2011.txt"};
  const double m_densityUnits = 1. / SI::eV / SI::m3;
  const double m_IgammaUnits = 1. / utils::pow<2>(SI::eV) / utils::pow<2>(SI::m3);

 public:
  Dominguez2011PhotonField();
  double getPhotonDensity(double ePhoton, double z = 0.) const override;
  double I_gamma(double ePhoton, double z = 0.) const override;
};

}  // namespace photonfield
}  // namespace simprop

#endif