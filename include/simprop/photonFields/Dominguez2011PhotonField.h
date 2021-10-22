#ifndef SIMPROP_PHOTONFIELDS_DOMINGUEZ2011_H
#define SIMPROP_PHOTONFIELDS_DOMINGUEZ2011_H

#include <string>
#include <vector>

#include "simprop/photonFields/AbstractPhotonField.h"

namespace simprop {
namespace photonfield {

class Dominguez2011Field : public AbstractField {
 public:
  Dominguez2011Field();
  double getPhotonDensity(double ePhoton, double z = 0.) const;
};

}  // namespace photonfield
}  // namespace simprop

#endif