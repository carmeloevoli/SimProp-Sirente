#include "simprop/pid.h"

#include <map>

namespace simprop {

PID getPidNucleus(int Z, unsigned int A) {
  assert(A >= 0);
  return PID(1000000000 + 10 * Z + 10000 * A);
}

const std::map<PID, std::string> pidNames = {{photon, "photon"},
                                             {neutrino_e, "nu_e"},
                                             {antineutrino_e, "antinu_e"},
                                             {neutrino_mu, "nu_mu"},
                                             {antineutrino_mu, "antinu_mu"},
                                             {electron, "electron"},
                                             {positron, "positron"},
                                             {proton, "proton"},
                                             {neutron, "neutron"},
                                             {antiproton, "antiproton"},
                                             {He4, "He4"},
                                             {C12, "C12"},
                                             {N14, "N14"},
                                             {O16, "O16"},
                                             {Fe56, "Fe56"}};

std::string getPidNames(const PID& pid) {
  auto it = pidNames.find(pid);
  if (it != pidNames.end())
    return it->second;
  else
    throw std::invalid_argument("pid name not found");
}

}  // namespace simprop
