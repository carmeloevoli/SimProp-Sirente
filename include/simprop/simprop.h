#ifndef SIMPROP_SIMPROP_H
#define SIMPROP_SIMPROP_H

#include <iostream>
#include <memory>
#include <vector>

#include "simprop/params.h"
#include "simprop/particle.h"
#include "simprop/photonFields/AbstractPhotonField.h"
#include "simprop/utils/random.h"

using RandomNumberGenerator = simprop::utils::RNG<double>;

namespace simprop {

class SimProp {
 private:
  const Params& m_params;
  size_t m_size;
  RandomNumberGenerator m_rng = utils::RNG<double>(1234);
  ParticleStack m_particles;
  std::vector<std::shared_ptr<photonfield::AbstractField> > m_photonFields;

 private:
  void printStateRanges() const;

 public:
  explicit SimProp(const Params& params);
  virtual ~SimProp();

  void buildInitialStates();
  void dumpParticles(std::string filename) const;

  void buildPhotonFields();
  void dumpPhotonFields() const;

  void run();
};

}  // namespace simprop

#endif  // SIMPROP_SIMPROP_H
