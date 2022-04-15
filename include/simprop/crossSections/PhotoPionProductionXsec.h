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
  size_t m_sSize = 9999;
  std::vector<double> m_sEnergies;
  std::vector<double> m_sigma;
  std::vector<double> m_phi;

 public:
  PhotoPionProductionXsec();
  virtual ~PhotoPionProductionXsec() = default;
  double getAtEpsPrime(double epsPrime) const override;
  double getAtS(double s) const override;
  double getPhiAtS(double s) const override;
  double getPhotonEnergyThreshold() const override;

 protected:
  void loadDataFile();
};

}  // namespace xsecs
}  // namespace simprop

#endif