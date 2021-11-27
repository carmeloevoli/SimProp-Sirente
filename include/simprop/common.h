#ifndef SIMPROP_COMMON_H
#define SIMPROP_COMMON_H

#include <utility>

#include "simprop/units.h"

namespace simprop {

double GetRndEnergy(std::pair<double, double> energyRange, double slope, double r);
double GetRndRedshift(std::pair<double, double> redshiftRange, int evolutionIndex, double r);

inline double energyToFrequency(double energy) { return energy / SI::hPlanck; }
inline double energyToWavelenght(double energy) { return SI::hPlanck * SI::cLight / energy; }

}  // namespace simprop

#endif