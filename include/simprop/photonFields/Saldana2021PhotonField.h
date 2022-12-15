#ifndef SIMPROP_PHOTONFIELDS_SALDANA2021_H
#define SIMPROP_PHOTONFIELDS_SALDANA2021_H

#include <string>
#include <vector>

#include "simprop/core/units.h"
#include "simprop/photonFields/LookupTablePhotonField.h"

namespace simprop {
namespace photonfields {

class Saldana2021PhotonField final : public LookupTablePhotonField {
 public:
  Saldana2021PhotonField() : LookupTablePhotonField(37, 250, "EBL_Saldana2021.txt") {}
};

class Saldana2021LowerPhotonField final : public LookupTablePhotonField {
 public:
  Saldana2021LowerPhotonField() : LookupTablePhotonField(37, 250, "EBL_lower_Saldana2021.txt") {}
};

class Saldana2021UpperPhotonField final : public LookupTablePhotonField {
 public:
  Saldana2021UpperPhotonField() : LookupTablePhotonField(37, 250, "EBL_upper_Saldana2021.txt") {}
};

}  // namespace photonfields
}  // namespace simprop

#endif