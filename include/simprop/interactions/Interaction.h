#ifndef SIMPROP_INTERACTIONS_INTERACTION_H
#define SIMPROP_INTERACTIONS_INTERACTION_H

#include <memory>

#include "simprop/crossSections/CrossSection.h"
#include "simprop/particle.h"

namespace simprop {
namespace interactions {

class Interaction {
 protected:
  std::shared_ptr<xsecs::CrossSection> m_sigma;

 public:
  Interaction(const std::shared_ptr<xsecs::CrossSection>& sigma) : m_sigma(sigma) {}
  virtual ~Interaction() = default;
  virtual double rate(PID pid, double Gamma, double z = 0) const = 0;
  virtual std::vector<Particle> finalState(const Particle& particle) const = 0;
};

}  // namespace interactions
}  // namespace simprop

#endif
