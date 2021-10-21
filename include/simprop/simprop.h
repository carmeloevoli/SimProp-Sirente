#ifndef SIMPROP_SIMPROP_H
#define SIMPROP_SIMPROP_H

#include <iostream>
#include <memory>
#include <vector>

#include "simprop/Units.h"
#include "simprop/params.h"
#include "simprop/photonFields/photonField.h"
#include "simprop/pid.h"
#include "simprop/utils/random.h"

using RandomNumberGenerator = simprop::utils::RNG<double>;

namespace simprop {

struct Particle {
  PID pid;
  double z;
  double E;

  friend std::ostream& operator<<(std::ostream& os, const Particle& p) {
    auto n = getPidNames(p.pid);
    return os << n << " " << std::setw(9) << p.z << " " << std::setw(16) << p.E / SI::eV;
  }
};

class SimProp {
 private:
  const Params& m_params;
  size_t m_size;
  RandomNumberGenerator m_rng = utils::RNG<double>(1234);
  std::vector<Particle> m_primaries;
  std::vector<photonfield::AbstractField> m_photonFields;

 private:
  void printRanges() const;

 public:
  explicit SimProp(const Params& params);
  virtual ~SimProp();

  void buildInitialStates();
  void dumpPrimaryParticles(std::string filename);

  void buildPhotonFields();
  void dumpPhotonFields(std::string filename);

  void run();
};

}  // namespace simprop

#endif  // SIMPROP_SIMPROP_H
