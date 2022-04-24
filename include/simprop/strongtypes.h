#ifndef SIMPROP_STRONGTYPES_H
#define SIMPROP_STRONGTYPES_H

#include "NamedType/named_type.hpp"

namespace simprop {

using RandomNumber = fluent::NamedType<double, struct RandomNumberTag, fluent::Comparable>;
using Seed = fluent::NamedType<long unsigned int, struct SeedTag>;
using Redshift = fluent::NamedType<double, struct RedshiftTag>;
using LorentzFactor = fluent::NamedType<double, struct LorentzFactorTag>;

}  // namespace simprop

#endif