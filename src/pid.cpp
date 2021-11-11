#include "simprop/pid.h"

#include <map>

namespace simprop {

PID getPidNucleus(const int& Z, const int& A) {
  if (A < 0 || Z > A) throw std::invalid_argument("invalid arguments for nucleus PID");
  return PID(1000000000 + 10 * Z + 10000 * A);
}

bool isNucleus(const PID& pid) { return (pid >= PID(1000009990)); }

int getNucleusMassNumber(const PID& pid) {
  if (!isNucleus(pid)) throw std::invalid_argument(getPidNames(pid) + " is not a nucleus");
  if (pid == neutron || pid == antiproton)
    return 1;
  else
    return ((int)pid / 10000) % 1000;
}

int getNucleusChargeNumber(const PID& pid) {
  if (!isNucleus(pid)) throw std::invalid_argument(getPidNames(pid) + " is not a nucleus");
  if (pid == neutron)
    return 0;
  else if (pid == antiproton)
    return -1;
  else
    return ((int)pid / 10) % 1000;
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
