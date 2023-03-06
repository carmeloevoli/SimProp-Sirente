#ifndef SIMPROP_XSECS_PHOTOPIONXSECS_H
#define SIMPROP_XSECS_PHOTOPIONXSECS_H

#include <string>

#include "simprop/crossSections/CrossSection.h"
#include "simprop/utils/lookupContainers.h"

namespace simprop {
namespace xsecs {

class PhotoPionXsec final : public CrossSection {
 private:
  const std::string m_filename = "data/xsecs_photopion_proton_sophia.txt";
  utils::LookupArray<2849> m_proton_sigma;
  utils::LookupArray<2849> m_proton_phi;
  utils::LookupArray<2850> m_neutron_sigma;
  utils::LookupArray<2850> m_neutron_phi;

 public:
  PhotoPionXsec();
  virtual ~PhotoPionXsec() = default;
  double getAtEpsPrime(PID pid, double eps) const override;
  double getPhotonEnergyThreshold() const override;

  double getAtS(PID pid, double s) const;
  double getPhiAtS(PID pid, double s) const;

 private:
  double getProtonXsec(double s) const;
  double getNeutronXsec(double s) const;
};

}  // namespace xsecs
}  // namespace simprop

#endif