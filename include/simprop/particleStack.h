#ifndef SIMPROP_PARTICLESTACK_H
#define SIMPROP_PARTICLESTACK_H

#include "simprop/particle.h"
#include "simprop/utils/random.h"

using RandomNumberGenerator = simprop::utils::RNG<double>;

namespace simprop {

class ParticleStack {
 protected:
  std::vector<Particle> m_stack;
  size_t m_size;
  RandomNumberGenerator m_rng = utils::RNG<double>(1234);

 public:
  ParticleStack() {}
  virtual ~ParticleStack() = default;

  void buildInitialStates();
  void printStateRanges() const;
};

}  // namespace simprop

#endif  // SIMPROP_PARTICLESTACK_H