#include "simprop/pid.h"

#include <map>
#include <string>

namespace simprop {

PID getPidNucleus(const int& Z, const int& A) {
  if (A < 0 || Z > A) throw std::invalid_argument("invalid arguments for nucleus PID");
  return PID(1000000000 + 10 * Z + 10000 * A);
}

bool isNucleus(const PID& pid) { return (pid >= PID(1000009990)); }

int getNucleusMassNumber(const PID& pid) {
  if (!isNucleus(pid)) throw std::invalid_argument(getPidName(pid) + " is not a nucleus");
  if (pid == neutron || pid == antiproton)
    return 1;
  else
    return ((int)pid / 10000) % 1000;
}

int getNucleusChargeNumber(const PID& pid) {
  if (!isNucleus(pid)) throw std::invalid_argument(getPidName(pid) + " is not a nucleus");
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
                                             {deuterium, "D"},
                                             {He3, "He3"},
                                             {He4, "He4"},
                                             {Be9, "Be9"},
                                             {B10, "B10"},
                                             {B11, "B11"},
                                             {C12, "C12"},
                                             {C13, "C13"},
                                             {N14, "N14"},
                                             {N15, "N15"},
                                             {O16, "O16"},
                                             {O17, "O17"},
                                             {O18, "O18"},
                                             {F19, "F19"},
                                             {Ne20, "Ne20"},
                                             {Ne21, "Ne21"},
                                             {Ne22, "Ne22"},
                                             {Na23, "Na23"},
                                             {Mg24, "Mg24"},
                                             {Mg25, "Mg25"},
                                             {Mg26, "Mg26"},
                                             {Al27, "Al27"},
                                             {Si28, "Si28"},
                                             {Si29, "Si29"},
                                             {Si30, "Si30"},
                                             {P31, "P31"},
                                             {S32, "S32"},
                                             {S33, "S33"},
                                             {S34, "S34"},
                                             {Cl35, "Cl35"},
                                             {Ar36, "Ar36"},
                                             {Cl37, "Cl37"},
                                             {Ar38, "Ar38"},
                                             {K39, "K39"},
                                             {Ca40, "Ca40"},
                                             {K41, "K41"},
                                             {Ca42, "Ca42"},
                                             {Ca43, "Ca43"},
                                             {Ca44, "Ca44"},
                                             {Sc45, "Sc45"},
                                             {Ti46, "Ti46"},
                                             {Ti47, "Ti47"},
                                             {Ti48, "Ti48"},
                                             {Ti49, "Ti49"},
                                             {Cr50, "Cr50"},
                                             {V51, "V51"},
                                             {Cr52, "Cr52"},
                                             {Cr53, "Cr53"},
                                             {Fe54, "Fe54"},
                                             {Mn55, "Mn55"},
                                             {Fe56, "Fe56"}};

std::string getPidName(const PID& pid) {
  auto it = pidNames.find(pid);
  if (it != pidNames.end())
    return it->second;
  else
    throw std::invalid_argument("pid name not found");
}

}  // namespace simprop
