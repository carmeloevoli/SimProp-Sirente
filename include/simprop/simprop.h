#ifndef SIMPROP_SIMPROP_H
#define SIMPROP_SIMPROP_H

#include <iostream>
#include <vector>

#include "simprop/Units.h"
#include "simprop/params.h"
#include "simprop/pid.h"
#include "simprop/utils/random.h"

using RandomNumberGenerator = simprop::utils::RNG<double>;

namespace simprop {

struct Particle {
  PID pid;
  double z;
  double E;

  friend std::ostream& operator<<(std::ostream& os, const Particle& p) {
    return os << getPidNames(p.pid) << " " << std::setw(9) << p.z << " " << std::setw(16)
              << p.E / SI::eV;
  }
};

class SimProp {
 private:
  const Params& m_params;
  size_t m_size;
  RandomNumberGenerator m_rng = utils::RNG<double>(1234);
  std::vector<Particle> m_primaries;

 private:
  void printRanges() const;

 public:
  explicit SimProp(const Params& params);
  virtual ~SimProp();

  void buildInitialStates();
  void dumpPrimaryParticles(std::string filename);
  void run();
};

}  // namespace simprop

#endif  // SIMPROP_SIMPROP_H