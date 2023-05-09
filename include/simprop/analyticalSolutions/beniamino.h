#ifndef SIMPROP_ANALYTICALSOLUTIONS_BENIAMINO_H
#define SIMPROP_ANALYTICALSOLUTIONS_BENIAMINO_H

#include "simprop/core/cosmology.h"
#include "simprop/energyLosses/ContinuousLosses.h"
#include "simprop/utils/lookupContainers.h"

namespace simprop {
namespace solutions {

struct SourceParams {
  double injSlope;
  double evolutionIndex;
  double expCutoff;
  double zMax;
};

class Beniamino {
 public:
  Beniamino(const SourceParams& params, const std::shared_ptr<cosmo::Cosmology>& cosmology,
            const std::vector<std::shared_ptr<losses::ContinuousLosses>>& losses);
  Beniamino& doCaching();

  virtual ~Beniamino() = default;

  double generationEnergy(double E, double zInit, double zFinal, double relError = 1e-3) const;
  double dilationFactor(double E, double zInit, double zFinal, double relError = 1e-3) const;
  double computeFlux(double E, double zObs, double relError = 1e-2) const;
  // double computeFluxUnm(double E, double zMax, double relError = 1e-3) const;

 public:
  double dbdE(double E) const;
  double beta(double E) const;
  double getMaxRedshift() const { return m_zMax; }
  const std::shared_ptr<cosmo::Cosmology>& getCosmology() const { return m_cosmology; }

 protected:
  std::shared_ptr<cosmo::Cosmology> m_cosmology;
  std::vector<std::shared_ptr<losses::ContinuousLosses>> m_losses;
  utils::LookupArray<10000> m_lossesLookup;

  const double m_sourceEmissivity{1e45 * SI::erg / SI::Mpc3 / SI::year};
  const double m_maxEnergy{1e23 * SI::eV};
  const double m_minEnergy{1e17 * SI::eV};

  double m_zMax{3.};
  double m_injSlope{2.6};
  double m_evolutionIndex{0.};
  double m_expCutoff{-1.};
  bool m_doCaching{false};
};

}  // namespace solutions
}  // namespace simprop

#endif