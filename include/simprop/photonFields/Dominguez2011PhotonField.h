#ifndef SIMPROP_PHOTONFIELDS_DOMINGUEZ2011_H
#define SIMPROP_PHOTONFIELDS_DOMINGUEZ2011_H

#include <string>
#include <vector>

#include "simprop/photonFields/LookupTablePhotonField.h"
#include "simprop/units.h"

namespace simprop {
namespace photonfields {

// enum class EblModel { MEAN, UPPER, LOWER };

class Dominguez2011PhotonField final : public LookupTablePhotonField {
 public:
  Dominguez2011PhotonField() : LookupTablePhotonField(18, 50, "EBL_Dominguez2011.txt") {}
};

class Dominguez2011LowerPhotonField final : public LookupTablePhotonField {
 public:
  Dominguez2011LowerPhotonField() : LookupTablePhotonField(18, 50, "EBL_lower_Dominguez2011.txt") {}
};

class Dominguez2011UpperPhotonField final : public LookupTablePhotonField {
 public:
  Dominguez2011UpperPhotonField() : LookupTablePhotonField(18, 50, "EBL_upper_Dominguez2011.txt") {}
};

}  // namespace photonfields
}  // namespace simprop

#endif