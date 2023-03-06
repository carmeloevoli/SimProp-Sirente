// Copyright 2023 SimProp-dev [MIT License]
#ifndef SIMPROP_CROSSSECTIONS_PHOTODISINTEGRATIONTALYSXSECS_H_
#define SIMPROP_CROSSSECTIONS_PHOTODISINTEGRATIONTALYSXSECS_H_

#include <string>

#include "simprop/crossSections/CrossSection.h"

namespace simprop {
namespace xsecs {

class PhotoDisintegrationTalysXsec final : public CrossSection {
 public:
  PhotoDisintegrationTalysXsec();
  virtual ~PhotoDisintegrationTalysXsec() = default;
  double getAtS(PID pid, double s) const override;
  double getPhiAtS(PID pid, double s) const override;
  double getPhotonEnergyThreshold() const override;
};

}  // namespace xsecs
}  // namespace simprop

#endif  // SIMPROP_CROSSSECTIONS_PHOTODISINTEGRATIONTALYSXSECS_H_
