#ifndef SIMPROP_PID_H
#define SIMPROP_PID_H

#include <cassert>
#include <stdexcept>

#include "NamedType/named_type.hpp"

namespace simprop {

using PID = fluent::NamedType<long int, struct PidTag, fluent::Comparable, fluent::Hashable>;

PID getPidNucleus(const int& Z, const int& A);
bool pidIsNucleus(const PID& pid);
bool pidIsNucleon(const PID& pid);
int getPidNucleusMassNumber(const PID& pid);
int getPidNucleusCharge(const PID& pid);
double getPidMass(const PID& pid);
std::string getPidName(const PID& pid);
PID removeNucleon(const PID& pid, const PID& nucleon);

static const PID photon = PID(0);
static const PID neutrino_e = PID(1);
static const PID antineutrino_e = PID(2);
static const PID neutrino_mu = PID(3);
static const PID antineutrino_mu = PID(4);
static const PID pionNeutral = PID(5);
static const PID pionPlus = PID(6);
static const PID pionMinus = PID(7);
static const PID electron = PID(8);
static const PID positron = PID(9);

static const PID neutron = getPidNucleus(0, 1);
static const PID proton = getPidNucleus(1, 1);
static const PID antiproton = getPidNucleus(-1, 1);
static const PID deuterium = getPidNucleus(1, 2);
static const PID He3 = getPidNucleus(2, 3);
static const PID He4 = getPidNucleus(2, 4);
static const PID Be9 = getPidNucleus(4, 9);
static const PID B10 = getPidNucleus(5, 10);
static const PID B11 = getPidNucleus(5, 11);
static const PID C12 = getPidNucleus(6, 12);
static const PID C13 = getPidNucleus(6, 13);
static const PID N14 = getPidNucleus(7, 14);
static const PID N15 = getPidNucleus(7, 15);
static const PID O16 = getPidNucleus(8, 16);
static const PID O17 = getPidNucleus(8, 17);
static const PID O18 = getPidNucleus(8, 18);
static const PID F19 = getPidNucleus(9, 19);
static const PID Ne20 = getPidNucleus(10, 20);
static const PID Ne21 = getPidNucleus(10, 21);
static const PID Ne22 = getPidNucleus(10, 22);
static const PID Na23 = getPidNucleus(11, 23);
static const PID Mg24 = getPidNucleus(12, 24);
static const PID Mg25 = getPidNucleus(12, 25);
static const PID Mg26 = getPidNucleus(12, 26);
static const PID Al27 = getPidNucleus(13, 27);
static const PID Si28 = getPidNucleus(14, 28);
static const PID Si29 = getPidNucleus(14, 29);
static const PID Si30 = getPidNucleus(14, 30);
static const PID P31 = getPidNucleus(15, 31);
static const PID S32 = getPidNucleus(16, 32);
static const PID S33 = getPidNucleus(16, 33);
static const PID S34 = getPidNucleus(16, 34);
static const PID Cl35 = getPidNucleus(17, 35);
static const PID Ar36 = getPidNucleus(18, 36);
static const PID Cl37 = getPidNucleus(17, 37);
static const PID Ar38 = getPidNucleus(18, 38);
static const PID K39 = getPidNucleus(19, 39);
static const PID Ca40 = getPidNucleus(20, 40);
static const PID K41 = getPidNucleus(19, 41);
static const PID Ca42 = getPidNucleus(20, 42);
static const PID Ca43 = getPidNucleus(20, 43);
static const PID Ca44 = getPidNucleus(20, 44);
static const PID Sc45 = getPidNucleus(21, 45);
static const PID Ti46 = getPidNucleus(22, 46);
static const PID Ti47 = getPidNucleus(22, 47);
static const PID Ti48 = getPidNucleus(22, 48);
static const PID Ti49 = getPidNucleus(22, 49);
static const PID Cr50 = getPidNucleus(24, 50);
static const PID V51 = getPidNucleus(23, 51);
static const PID Cr52 = getPidNucleus(24, 52);
static const PID Cr53 = getPidNucleus(24, 53);
static const PID Fe54 = getPidNucleus(26, 54);
static const PID Mn55 = getPidNucleus(25, 55);
static const PID Fe56 = getPidNucleus(26, 56);

}  // namespace simprop

#endif  // SIMPROP_PID_H