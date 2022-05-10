#ifndef SIMPROP_XSECS_PHOTOPIONPRODUCTION_H
#define SIMPROP_XSECS_PHOTOPIONPRODUCTION_H

#include <string>

#include "simprop/crossSections/CrossSection.h"
#include "simprop/utils/lookupContainers.h"

namespace simprop {
namespace xsecs {

class PhotoPionProductionXsec final : public CrossSection {
 protected:
  const std::string m_filename = "data/xsec_ppp.txt";
  utils::LookupArray<9999> m_sigma;
  utils::LookupArray<9999> m_phi;

 public:
  PhotoPionProductionXsec();
  virtual ~PhotoPionProductionXsec() = default;
  double getAtEpsPrime(double epsPrime) const override;
  double getAtS(double s) const override;
  double getPhiAtS(double s) const override;
  double getPhotonEnergyThreshold() const override;
};

}  // namespace xsecs
}  // namespace simprop

#endif