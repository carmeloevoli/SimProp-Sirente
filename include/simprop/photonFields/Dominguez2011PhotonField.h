// Copyright 2023 SimProp-dev [MIT License]
#ifndef SIMPROP_PHOTONFIELDS_DOMINGUEZ2011PHOTONFIELD_H_
#define SIMPROP_PHOTONFIELDS_DOMINGUEZ2011PHOTONFIELD_H_

#include <string>
#include <vector>

#include "simprop/core/units.h"
#include "simprop/photonFields/LookupTablePhotonField.h"

namespace simprop {
namespace photonfields {

class Dominguez2011PhotonField final : public LookupTablePhotonField {
 public:
  Dominguez2011PhotonField() : LookupTablePhotonField(18, 50, "ebl_Dominguez2011_fiducial.txt") {}
};

class Dominguez2011LowerPhotonField final : public LookupTablePhotonField {
 public:
  Dominguez2011LowerPhotonField() : LookupTablePhotonField(18, 50, "ebl_Dominguez2011_lower.txt") {}
};

class Dominguez2011UpperPhotonField final : public LookupTablePhotonField {
 public:
  Dominguez2011UpperPhotonField() : LookupTablePhotonField(18, 50, "ebl_Dominguez2011_upper.txt") {}
};

}  // namespace photonfields
}  // namespace simprop

#endif  // SIMPROP_PHOTONFIELDS_DOMINGUEZ2011PHOTONFIELD_H_
