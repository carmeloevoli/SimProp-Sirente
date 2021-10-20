#ifndef SIMPROP_PID_H
#define SIMPROP_PID_H

#include <cassert>
#include <stdexcept>

namespace simprop {

enum class PID : long int {};

PID getPidNucleus(int Z, unsigned int A);
std::string getPidNames(const PID& pid);

static const PID photon = PID(0);
static const PID neutrino_e = PID(1);
static const PID antineutrino_e = PID(2);
static const PID neutrino_mu = PID(3);
static const PID antineutrino_mu = PID(4);
static const PID electron = PID(5);
static const PID positron = PID(6);

static const PID neutron = getPidNucleus(0, 1);
static const PID proton = getPidNucleus(1, 1);
static const PID antiproton = getPidNucleus(-1, 1);
static const PID He4 = getPidNucleus(2, 4);
static const PID C12 = getPidNucleus(6, 12);
static const PID N14 = getPidNucleus(7, 14);
static const PID O16 = getPidNucleus(8, 16);
static const PID Fe56 = getPidNucleus(26, 56);

}  // namespace simprop

#endif  // SIMPROP_PID_H