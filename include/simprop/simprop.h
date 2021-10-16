#ifndef SIMPROP_SIMPROP_H
#define SIMPROP_SIMPROP_H

#include <vector>

#include "simprop/Units.h"
#include "simprop/params.h"
#include "simprop/pid.h"
#include "simprop/random.h"

using RandomNumberGenerator = simprop::utils::RNG<double>;

namespace simprop {

class PrimaryParticle {
 private:
  PID m_pid;
  double m_z;
  double m_E;

 public:
  explicit PrimaryParticle(PID pid, double z, double E) : m_pid(pid), m_z(z), m_E(E) {}
  virtual ~PrimaryParticle(){};
};

class SimProp {
 private:
  const Params& m_params;
  size_t m_size;
  RandomNumberGenerator m_rng = utils::RNG<double>(12345);
  std::vector<PrimaryParticle> m_primaries;

 public:
  explicit SimProp(const Params& params);
  virtual ~SimProp();

  void buildInitialStates();
  void run();
};

}  // namespace simprop

#endif  // SIMPROP_SIMPROP_H
