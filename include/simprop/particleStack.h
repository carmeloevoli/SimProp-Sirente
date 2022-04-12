#ifndef SIMPROP_PARTICLESTACK_H
#define SIMPROP_PARTICLESTACK_H

#include <memory>

#include "simprop/particle.h"
#include "simprop/utils/random.h"

using Range = std::pair<double, double>;

namespace simprop {

class ParticleStack {
 public:
  using stack = std::vector<Particle>;
  using iterator = typename stack::iterator;
  using const_iterator = typename stack::const_iterator;

 protected:
  PID m_pid;
  size_t m_size;
  stack m_particles;
  RandomNumberGenerator& m_rng;

 public:
  iterator begin() { return m_particles.begin(); }
  iterator end() { return m_particles.end(); }
  const_iterator begin() const { return m_particles.begin(); }
  const_iterator end() const { return m_particles.end(); }

 public:
  explicit ParticleStack(PID pid, int nParticles, RandomNumberGenerator& rng);
  virtual ~ParticleStack() = default;

  void buildInitialState(Range zRange, Range gammaRange, double slope);
  void buildSingleParticleStack(double z, double Gamma);
  void buildMultipleParticleStack(Range zRange, Range gammaRange, double slope);

  Range getRedshiftRange() const;
  Range getGammaRange() const;
};

}  // namespace simprop

#endif  // SIMPROP_PARTICLESTACK_H