#ifndef SIMPROP_ANALYTICALSOLUTIONS_COSMONUS_H
#define SIMPROP_ANALYTICALSOLUTIONS_COSMONUS_H

#include <memory>

#include "simprop/core/cosmology.h"
#include "simprop/crossSections/KelnerAharonian2008.h"
#include "simprop/photonFields/CmbPhotonField.h"
#include "simprop/photonFields/Saldana2021PhotonField.h"
#include "simprop/utils/lookupContainers.h"

namespace simprop {
namespace solutions {

class CosmoNeutrinos {
 protected:
  std::shared_ptr<cosmo::Cosmology> m_cosmology;
  std::shared_ptr<photonfields::PhotonField> m_cmb;
  std::shared_ptr<KelnerAharonian2008::NeutrinoProductionSpectrum> m_nuSpec;

  utils::LookupTable<50, 11> m_Jp;

 public:
  CosmoNeutrinos();
  virtual ~CosmoNeutrinos() = default;

  double computeNeutrinoFlux(double Enu, double zMax, size_t N = 8) const;
  double getProtonFlux(double Ep, double z) const;

 protected:
  double I_deps(double EnuObs, double Ep, double z, size_t N = 10) const;
  double I_dEp(double EnuObs, double z, size_t N = 10) const;
};

}  // namespace solutions
}  // namespace simprop

#endif  // SIMPROP_ANALYTICALSOLUTIONS_COSMONUS_H