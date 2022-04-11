#include "simprop/pid.h"

#include <map>
#include <string>

#include "simprop/units.h"

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

int getNucleusCharge(const PID& pid) {
  if (!isNucleus(pid)) throw std::invalid_argument(getPidName(pid) + " is not a nucleus");
  if (pid == neutron)
    return 0;
  else if (pid == antiproton)
    return -1;
  else
    return ((int)pid / 10) % 1000;
}

static const std::map<PID, std::string> pidNames = {
    {photon, "photon"},       {neutrino_e, "nu_e"},           {antineutrino_e, "antinu_e"},
    {neutrino_mu, "nu_mu"},   {antineutrino_mu, "antinu_mu"}, {electron, "electron"},
    {positron, "positron"},   {pionNeutral, "pion_0"},        {pionPlus, "pion_plus"},
    {pionMinus, "pion_minus"}};

static const std::map<int, std::string> chargeToName = {
    {1, "H"},   {2, "He"}, {3, "Li"},  {4, "Be"},  {5, "B"},   {6, "C"},   {7, "N"},
    {8, "O"},   {9, "F"},  {10, "Ne"}, {11, "Na"}, {12, "Mg"}, {13, "Al"}, {14, "Si"},
    {15, "P"},  {16, "S"}, {17, "Cl"}, {18, "Ar"}, {19, "K"},  {20, "Ca"}, {21, "Sc"},
    {22, "Ti"}, {23, "V"}, {24, "Cr"}, {25, "Mn"}, {26, "Fe"}, {27, "Co"}, {28, "Ni"}};

std::string getNucleusName(const PID& pid) {
  auto A = getNucleusMassNumber(pid);
  auto Z = getNucleusCharge(pid);
  auto it = chargeToName.find(Z);
  if (it != chargeToName.end())
    return it->second + std::to_string(A);
  else
    throw std::invalid_argument("pid name not found");
}

double getMassFromPid(const PID& pid) {
  if (isNucleus(pid)) {
    auto A = (double)getNucleusMassNumber(pid);
    auto Z = (double)getNucleusCharge(pid);
    return (A - Z) * SI::neutronMassC2 + Z * SI::protonMassC2;
  } else if (pid == positron || pid == electron)
    return SI::electronMassC2;
  else if (pid == pionNeutral || pid == pionMinus || pid == pionPlus)
    return SI::pionMassC2;
  else
    throw std::invalid_argument("mass not available for this pid");
}

std::string getPidName(const PID& pid) {
  if (pid == proton) return "proton";
  if (pid == neutron) return "neutron";
  if (pid == antiproton) return "antiproton";
  if (pid == deuterium) return "deuterium";
  if (isNucleus(pid)) return getNucleusName(pid);
  auto it = pidNames.find(pid);
  if (it != pidNames.end())
    return it->second;
  else
    throw std::invalid_argument("pid name not found");
}

}  // namespace simprop
