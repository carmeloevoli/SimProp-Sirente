#ifndef SIMPROP_PHOTONFIELDS_GILMORE2012_H
#define SIMPROP_PHOTONFIELDS_GILMORE2012_H

#include <string>
#include <vector>

#include "simprop/core/units.h"
#include "simprop/photonFields/LookupTablePhotonField.h"

namespace simprop {
namespace photonfields {

class Gilmore2012PhotonField final : public LookupTablePhotonField {
 public:
  Gilmore2012PhotonField() : LookupTablePhotonField(20, 101, "ebl_Gilmore2012_fiducial.txt") {}
};

}  // namespace photonfields
}  // namespace simprop

#endif