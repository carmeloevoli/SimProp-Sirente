#ifndef SIMPROP_PARTICLE_H
#define SIMPROP_PARTICLE_H

#include "simprop/pid.h"
#include "simprop/units.h"

namespace simprop {

class Particle {
  struct State {
    double z;
    double E;
  };

 protected:
  PID m_pid;
  State m_origin;
  State m_now;

 public:
  Particle(PID pid, double z, double E) : m_pid(pid), m_origin({z, E}), m_now({z, E}) {}

  const State getNow() const { return m_now; }
  State& getNow() { return m_now; }

  const double getRedshift() const { return m_now.z; }
  const double getEnergy() const { return m_now.E; }
  const PID getPid() const { return m_pid; }

  friend std::ostream& operator<<(std::ostream& os, const Particle& p) {
    auto n = getPidName(p.m_pid);
    auto z = p.m_now.z;
    auto E = p.m_now.E / SI::eV;
    return os << n << "\t" << std::scientific << z << "\t" << E;
  }
};

}  // namespace simprop

#endif