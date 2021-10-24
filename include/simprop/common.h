#ifndef SIMPROP_COMMON_H
#define SIMPROP_COMMON_H

#include <utility>

#include "simprop/Units.h"

namespace simprop {

double GetRndEnergy(std::pair<double, double> energyRange, double r);
double GetRndRedshift(double maxRedshift, int evolutionIndex, double r);

inline double energyToFrequency(double energy) { return energy / SI::hPlanck; }

}  // namespace simprop

#endif