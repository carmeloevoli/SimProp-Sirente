#ifndef SIMPROP_XSECS_PHOTOPIONXSECS_H
#define SIMPROP_XSECS_PHOTOPIONXSECS_H

#include <string>

#include "simprop/crossSections/CrossSection.h"
#include "simprop/utils/lookupContainers.h"

namespace simprop {
namespace xsecs {

class PhotoPionProtonXsec final : public CrossSection {
 private:
  const std::string m_filename = "data/xsecs_photopion_proton_sophia.txt";
  utils::LookupArray<2849> m_sigma;
  utils::LookupArray<2849> m_phi;

 public:
  PhotoPionProtonXsec();

 public:
  virtual ~PhotoPionProtonXsec() = default;
  double getAtS(double s) const override;
  double getPhiAtS(double s) const override;
  double getPhotonEnergyThreshold() const override;
};

class PhotoPionNeutronXsec final : public CrossSection {
 private:
  const std::string m_filename = "data/xsecs_photopion_neutron_sophia.txt";
  utils::LookupArray<2850> m_sigma;
  utils::LookupArray<2850> m_phi;

 public:
  PhotoPionNeutronXsec();
  virtual ~PhotoPionNeutronXsec() = default;
  double getAtS(double s) const override;
  double getPhiAtS(double s) const override;
  double getPhotonEnergyThreshold() const override;
};

}  // namespace xsecs
}  // namespace simprop

#endif