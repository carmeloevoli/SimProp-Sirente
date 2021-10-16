#ifndef SIMPROP_COMMON_H
#define SIMPROP_COMMON_H

#include <utility>

#include "simprop/Units.h"

namespace simprop {

double GetRndEnergy(std::pair<double, double> energyRange, double r);
double GetRndRedshift(double maxRedshift, int evolutionIndex, double r);

}  // namespace simprop

#endif