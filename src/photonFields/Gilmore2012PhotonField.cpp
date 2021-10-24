#include "simprop/photonFields/Gilmore2012PhotonField.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include "simprop/utils/misc.h"

namespace simprop {
namespace photonfield {

Gilmore2012Field::Gilmore2012Field() { m_fieldName = "Gilmore2012"; }

double Gilmore2012Field::getPhotonDensity(double ePhoton, double z) const {
  return table.get(ePhoton, z);
};

}  // namespace photonfield
}  // namespace simprop