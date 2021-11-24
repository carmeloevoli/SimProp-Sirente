#ifndef SIMPROP_PHOTONFIELDS_DOMINGUEZ2011_H
#define SIMPROP_PHOTONFIELDS_DOMINGUEZ2011_H

#include <string>
#include <vector>

#include "simprop/photonFields/AbstractPhotonField.h"
#include "simprop/utils/lookupTable.h"

namespace simprop {
namespace photonfield {

class Dominguez2011PhotonField : public AbstractPhotonField {
 protected:
  utils::LookupTable<18, 50> m_field{"data/EBL_Dominguez2011.txt"};

 public:
  Dominguez2011PhotonField();
  double getPhotonDensity(double ePhoton, double z = 0.) const override;
  double I_gamma(double ePhoton, double z = 0.) const override;
};

}  // namespace photonfield
}  // namespace simprop

#endif