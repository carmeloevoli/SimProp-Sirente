#ifndef SIMPROP_SIMPROP_H
#define SIMPROP_SIMPROP_H

#include <vector>

#include "simprop/Units.h"
#include "simprop/params.h"
#include "simprop/random.h"

using RandomNumberGenerator = simprop::utils::RNG<double>;
using InitialStates = std::pair<double, double>;
using PID = double;

namespace simprop {

class SimProp {
 private:
  const Params& m_params;
  RandomNumberGenerator m_rng = utils::RNG<double>(12345);
  std::vector<InitialStates> m_istates;
  PID pid;

 public:
  explicit SimProp(const Params& params);
  virtual ~SimProp();

  void buildInitialStates();
  void run();
};

}  // namespace simprop

#endif  // SIMPROP_SIMPROP_H
