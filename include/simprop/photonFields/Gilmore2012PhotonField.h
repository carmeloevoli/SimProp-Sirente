// Copyright 2023 SimProp-dev [MIT License]
#ifndef SIMPROP_PHOTONFIELDS_GILMORE2012PHOTONFIELD_H_
#define SIMPROP_PHOTONFIELDS_GILMORE2012PHOTONFIELD_H_

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

#endif  // SIMPROP_PHOTONFIELDS_GILMORE2012PHOTONFIELD_H_
