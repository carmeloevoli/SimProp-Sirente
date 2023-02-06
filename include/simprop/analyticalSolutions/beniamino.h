#ifndef SIMPROP_ANALYTICALSOLUTIONS_BENIAMINO_H
#define SIMPROP_ANALYTICALSOLUTIONS_BENIAMINO_H

#include "simprop/core/cosmology.h"
#include "simprop/crossSections/KelnerAharonian2008.h"
#include "simprop/energyLosses/BGG2006ContinuousLosses.h"
#include "simprop/energyLosses/PairProductionLosses.h"
#include "simprop/energyLosses/PhotoPionContinuousLosses.h"
#include "simprop/photonFields/CmbPhotonField.h"

namespace simprop {
namespace solutions {

class Beniamino {
 protected:
  std::shared_ptr<cosmo::Cosmology> m_cosmology;
  std::shared_ptr<photonfields::CMB> m_cmb;
  utils::LookupArray<4000> m_losses;

  const double m_sourceEmissivity = SI::GeV / SI::m3 / SI::sec;
  const double m_maxEnergy = 1e23 * SI::eV;
  double m_slope = 2.7;
  double m_sourceEvolution = 0.;
  double m_sourceCutoff = -1.;

  KelnerAharonian2008::NuMuSpectrum numu{};
  KelnerAharonian2008::AntiNuMuSpectrum antiNumu{};
  KelnerAharonian2008::NuElectronSpectrum nue{};
  KelnerAharonian2008::AntiNuElectronSpectrum antiNue{};

 public:
  Beniamino(bool doPhotoPion);
  virtual ~Beniamino() = default;

  double generationEnergy(double E, double zMax, double relError = 1e-3) const;
  double dilationFactor(double E, double zMax, double relError = 1e-3) const;
  double computeFlux(double E, double zObs, double zMax, double relError = 1e-3) const;
  double computeFluxUnm(double E, double zMax, double relError = 1e-3) const;
  double computeNeutrinoFlux(double Enu, double zMax, double relError = 1e-3) const;

 public:
  inline void setSlope(const double &slope) { m_slope = slope; };
  inline void setSourceEvolution(const double &m) { m_sourceEvolution = m; };
  inline void setSourceCutoff(const double &E) { m_sourceCutoff = E; };

 protected:
  double dbdE(double E) const;
  double findMaxRedshiftIntegral(double E, double zMax) const;
  double getNeutrinoSpectrum(double eta, double x) const;

 public:  // TODO make it private
  double I_dEpsilon(double EnuObserved, double Ep, double z) const;
  double I_dEp(double EnuObserved, double minEp, double maxEp, double z) const;
  // double I_z(double Enu, double zMax, double relError) const;
};

}  // namespace solutions
}  // namespace simprop

#endif