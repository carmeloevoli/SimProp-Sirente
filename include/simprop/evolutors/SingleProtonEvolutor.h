// Copyright 2023 SimProp-dev [MIT License]
#ifndef SIMPROP_EVOLUTORS_SINGLEPROTONEVOLUTORS_H
#define SIMPROP_EVOLUTORS_SINGLEPROTONEVOLUTORS_H

#include <memory>

#include "simprop/core/cosmology.h"
#include "simprop/energyLosses/ContinuousLosses.h"
#include "simprop/interactions/Interaction.h"
#include "simprop/particleStacks/Builder.h"
#include "simprop/utils/random.h"

namespace simprop {
namespace evolutors {

class SingleProtonEvolutor {
 public:
  SingleProtonEvolutor(RandomNumberGenerator& rng);
  virtual ~SingleProtonEvolutor() = default;

  void addCosmology(std::shared_ptr<cosmo::Cosmology> cosmology) { m_cosmology = cosmology; }
  void addLosses(std::vector<std::shared_ptr<losses::ContinuousLosses>> losses) {
    m_continuousLosses = losses;
  }
  void addInteractions(std::vector<std::shared_ptr<interactions::Interaction>> interactions) {
    m_interactions = interactions;
  }
  void run(ParticleStack& stack);

 protected:
  double computeDeltaGamma(const Particle& particle, double deltaRedshift) const;
  double computeLossesRedshiftInterval(const Particle& particle) const;
  double computeInteractionRedshiftInterval(const Particle& particle) const;
  double totalLosses(PID pid, double Gamma, double z) const;
  double totalRate(PID pid, double Gamma, double z) const;

 protected:
  const double deltaGammaCritical = 0.05;
  RandomNumberGenerator& m_rng;
  std::shared_ptr<cosmo::Cosmology> m_cosmology;
  std::vector<std::shared_ptr<losses::ContinuousLosses>> m_continuousLosses;
  std::vector<std::shared_ptr<interactions::Interaction>> m_interactions;
};

}  // namespace evolutors
}  // namespace simprop

#endif  // SIMPROP_EVOLUTORS_SINGLEPROTONEVOLUTORS_H
