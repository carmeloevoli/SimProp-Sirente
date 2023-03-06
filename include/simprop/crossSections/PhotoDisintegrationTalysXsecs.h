// Copyright 2023 SimProp-dev [MIT License]
#ifndef SIMPROP_CROSSSECTIONS_PHOTODISINTEGRATIONTALYSXSECS_H_
#define SIMPROP_CROSSSECTIONS_PHOTODISINTEGRATIONTALYSXSECS_H_

#include <string>

#include "simprop/crossSections/CrossSection.h"

namespace simprop {
namespace xsecs {

class PhotoDisintegrationTalysXsec final : public CrossSection {
 protected:
  std::string m_singleNucleonFilename = "data/";
  std::string m_alphaFilename = "data/";

 public:
  PhotoDisintegrationTalysXsec();
  virtual ~PhotoDisintegrationTalysXsec() = default;

  double getAtEpsPrime(PID pid, double eps) const override;
  double getPhotonEnergyThreshold() const override;
};

}  // namespace xsecs
}  // namespace simprop

#endif  // SIMPROP_CROSSSECTIONS_PHOTODISINTEGRATIONTALYSXSECS_H_
