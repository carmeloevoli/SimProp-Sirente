#ifndef SIMPROP_XSECS_PHOTOPIONPRODUCTION_H
#define SIMPROP_XSECS_PHOTOPIONPRODUCTION_H

#include <string>
#include <vector>

#include "simprop/crossSections/CrossSection.h"
#include "simprop/utils/lookupContainers.h"

namespace simprop {
namespace xsecs {

class PhotoPionProductionXsec final : public CrossSection {
 protected:
  const std::string m_filename = "data/xsec_ppp.txt";
  utils::LookupArray<1000> m_sigmas{m_filename};

 public:
  PhotoPionProductionXsec() {}
  virtual ~PhotoPionProductionXsec() = default;
  double getAtEpsPrime(PID pid, double epsPrime) const override;
  double getAtS(PID pid, double s) const override;
  double getPhotonEnergyThreshold() const override;
};

}  // namespace xsecs
}  // namespace simprop

#endif