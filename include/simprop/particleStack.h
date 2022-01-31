#ifndef SIMPROP_PARTICLESTACK_H
#define SIMPROP_PARTICLESTACK_H

#include "simprop/particle.h"
#include "simprop/utils/random.h"

using Range = std::pair<double, double>;

namespace simprop {

class ParticleStack {
  using stack = std::vector<Particle>;
  using iterator = typename stack::iterator;
  using const_iterator = typename stack::const_iterator;

 protected:
  PID m_pid;
  size_t m_size;
  RandomNumberGenerator m_rng = utils::RNG<double>(1234);
  stack m_particles;

 public:
  iterator begin() { return m_particles.begin(); }
  iterator end() { return m_particles.end(); }
  const_iterator begin() const { return m_particles.begin(); }
  const_iterator end() const { return m_particles.end(); }

 public:
  explicit ParticleStack(PID pid, int nParticles, int seed = 1234);
  virtual ~ParticleStack() = default;

  void buildInitialStates(Range zRange, Range eRange, double slope);
  Range getRedshiftRange() const;
  Range getEnergyRange() const;
};

}  // namespace simprop

#endif  // SIMPROP_PARTICLESTACK_H