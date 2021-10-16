#ifndef INCLUDE_SIMPROP_UNITS_H
#define INCLUDE_SIMPROP_UNITS_H

#include "simprop/utils.h"

using namespace simprop::utils;

namespace SI {
// from https://en.wikipedia.org/wiki/International_System_of_Units

// SI  base units
static const double meter = 1;
static const double second = 1;
static const double metre = 1;
static const double kilogram = 1;
static const double ampere = 1;
static const double kelvin = 1;

// SI derived units
static const double radian = meter / meter;
static const double hertz = 1 / second;
static const double joule = kilogram * pow<2>(meter / second);

// TODO

static const double eV = 1;
static const double GeV = 1e9;

}  // namespace SI

#endif  // INCLUDE_SIMPROP_UNITS_H