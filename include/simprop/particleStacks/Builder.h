#ifndef SIMPROP_PARTICLESTACK_H
#define SIMPROP_PARTICLESTACK_H

#include <memory>

#include "simprop/particle.h"
#include "simprop/utils/io.h"
#include "simprop/utils/random.h"

namespace simprop {

using ParticleStack = std::vector<Particle>;
using Range = std::pair<double, double>;

class Builder {
 protected:
  PID m_pid;
  size_t m_size;

 public:
  Builder(PID pid, size_t size = 1) : m_pid(pid), m_size(size) {}
  virtual ~Builder() = default;
  virtual ParticleStack build(RandomNumberGenerator& rng) const = 0;

  Range getRedshiftRange(const ParticleStack& stack) const;
  Range getGammaRange(const ParticleStack& stack) const;
};

}  // namespace simprop

#endif  // SIMPROP_PARTICLESTACK_H