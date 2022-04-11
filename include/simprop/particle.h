#ifndef SIMPROP_PARTICLE_H
#define SIMPROP_PARTICLE_H

#include <iomanip>

#include "simprop/pid.h"
#include "simprop/units.h"

namespace simprop {

class Particle {
  struct State {
    double z;
    double Gamma;
  };

 protected:
  PID m_pid;
  State m_origin;
  State m_now;

 public:
  Particle(PID pid, double z, double Gamma) : m_pid(pid), m_origin({z, Gamma}), m_now({z, Gamma}) {}

  const State getNow() const { return m_now; }
  State& getNow() { return m_now; }

  const double getRedshift() const { return m_now.z; }
  const double getGamma() const { return m_now.Gamma; }
  const double getEnergy() const { return m_now.Gamma * getMassFromPid(m_pid); }

  const PID getPid() const { return m_pid; }

  friend std::ostream& operator<<(std::ostream& os, const Particle& p) {
    auto pidName = getPidName(p.m_pid);
    auto z = p.m_now.z;
    auto Gamma = p.m_now.Gamma;
    return os << pidName << " " << std::fixed << z << " " << std::scientific << Gamma;
  }
};

}  // namespace simprop

#endif