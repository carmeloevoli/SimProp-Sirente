#ifndef SIMPROP_ANALYTICALSOLUTIONS_BZ_H
#define SIMPROP_ANALYTICALSOLUTIONS_BZ_H

#include <memory>

#include "simprop/core/cosmology.h"
#include "simprop/crossSections/KelnerAharonian2008.h"
#include "simprop/photonFields/CmbPhotonField.h"
#include "simprop/utils/lookupContainers.h"

namespace simprop {
namespace solutions {

class BzNeutrinos {
 protected:
  std::shared_ptr<cosmo::Cosmology> m_cosmology;
  std::shared_ptr<photonfields::CMB> m_cmb;
  std::shared_ptr<KelnerAharonian2008::NeutrinoProductionSpectrum> m_nuSpec;

  utils::LookupTable<100, 30> m_Jp;

 public:
  BzNeutrinos();
  virtual ~BzNeutrinos() = default;

  double computeNeutrinoFlux(double Enu, double zMax, double relError = 1e-3) const;
  double getProtonFlux(double Ep, double z) const;

 protected:
  double I_deps(double EnuObs, double Ep, double z) const;
  double I_dEp(double EnuObs, double z) const;
};

}  // namespace solutions
}  // namespace simprop

#endif  // SIMPROP_ANALYTICALSOLUTIONS_BZ_H