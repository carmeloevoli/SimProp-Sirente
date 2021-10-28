#ifndef SIMPROP_COSMOLOGY_COSMOLOGY_H
#define SIMPROP_COSMOLOGY_COSMOLOGY_H

namespace simprop {
namespace cosmo {

/* returns the hubble "constant" at z */
double H(double z);

/* returns hubble time, t_h = 1/H */
double t_H(double z);

/* returns the value of dt/dz at the redshift parameter z. */
double dtdz(double z);

/* returns the relative energy fraction lost by adiabatic losses */
double adiabaticLosses(double z_i, double z_f);

}  // namespace cosmo
}  // namespace simprop

#endif