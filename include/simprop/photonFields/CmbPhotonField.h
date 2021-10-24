#ifndef SIMPROP_PHOTONFIELDS_CMB_H
#define SIMPROP_PHOTONFIELDS_CMB_H

#include <string>

#include "simprop/photonFields/AbstractPhotonField.h"

namespace simprop {
namespace photonfield {

class BlackbodyPhotonField : public AbstractField {
 public:
  BlackbodyPhotonField(const std::string fieldName, const double blackbodyTemperature);
  double getPhotonDensity(double ePhoton, double z = 0.) const override;

 protected:
  double m_blackbodyTemperature;
};

class CMB : public BlackbodyPhotonField {
 public:
  CMB() : BlackbodyPhotonField("CMB", 2.73) {}
};

}  // namespace photonfield
}  // namespace simprop

#endif