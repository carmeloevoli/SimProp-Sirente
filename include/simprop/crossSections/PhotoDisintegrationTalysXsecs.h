// Copyright 2023 SimProp-dev [MIT License]
#ifndef SIMPROP_CROSSSECTIONS_PHOTODISINTEGRATIONTALYSXSECS_H_
#define SIMPROP_CROSSSECTIONS_PHOTODISINTEGRATIONTALYSXSECS_H_

#include <string>
#include <unordered_map>
#include <vector>

#include "simprop/crossSections/CrossSection.h"

namespace simprop {
namespace xsecs {

using XsecsMap = std::unordered_map<PID, std::vector<double>>;

class TalysChannel {
 protected:
  std::vector<double> m_energyAxis;
  XsecsMap m_xmap;

 public:
  TalysChannel();
  virtual ~TalysChannel() = default;
  void loadXsecMaps(const std::string filename);
  double get(PID pid, double eps) const;

 protected:
  void buildEnergyAxis();
};

class PhotoDisintegrationTalysXsec final : public CrossSection {
 protected:
  std::string m_singleNucleonFilename = "data/xsecs_photodisintegration_v2r4_singlenucleon.txt";
  std::string m_alphaFilename = "data/xsecs_photodisintegration_v2r4_alpha.txt";

  TalysChannel m_xsec_single;
  TalysChannel m_xsec_alpha;

 public:
  PhotoDisintegrationTalysXsec();
  virtual ~PhotoDisintegrationTalysXsec() = default;

  double getAtEpsPrime(PID pid, double eps) const override;
  double getPhotonEnergyThreshold() const override;
};

}  // namespace xsecs
}  // namespace simprop

#endif  // SIMPROP_CROSSSECTIONS_PHOTODISINTEGRATIONTALYSXSECS_H_
