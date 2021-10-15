#ifndef SIMPROP_SIMPROP_H
#define SIMPROP_SIMPROP_H

#include <vector>

#include "simprop/params.h"
#include "simprop/random.h"

namespace simprop {

using PID = double;

class SimProp {
 private:
  utils::RNG<double> m_rng = utils::RNG<double>(12345);
  std::vector<double> m_redshifts;
  std::vector<QEnergy> m_energies;
  PID pid;

 public:
  explicit SimProp(const Params& params);

  void init();
  void run();
};

}  // namespace simprop

#endif  // SIMPROP_SIMPROP_H
