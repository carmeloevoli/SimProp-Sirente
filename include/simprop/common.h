#ifndef SIMPROP_COMMON_H
#define SIMPROP_COMMON_H

#include <utility>

#include "simprop/Units.h"

namespace simprop {

double GetRndEnergy(std::pair<double, double> energyRange, double slope, double r);
double GetRndRedshift(std::pair<double, double> redshiftRange, int evolutionIndex, double r);

inline double energyToFrequency(double energy) { return energy / SI::hPlanck; }

}  // namespace simprop

#endif