#ifndef SIMPROP_PHOTONFIELDS_GILMORE2012_H
#define SIMPROP_PHOTONFIELDS_GILMORE2012_H

#include <string>
#include <vector>

#include "simprop/photonFields/AbstractPhotonField.h"
#include "simprop/photonFields/TabularPhotonField.h"

namespace simprop {
namespace photonfield {

class Gilmore2012Field : public AbstractField {
 public:
  Gilmore2012Field();
  double getPhotonDensity(double ePhoton, double z = 0.) const;

 protected:
  TabularPhotonField<100, 100> table{"ciccio"};
};

}  // namespace photonfield
}  // namespace simprop

#endif