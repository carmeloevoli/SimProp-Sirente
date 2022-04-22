#ifndef SIMPROP_PARTICLE_H
#define SIMPROP_PARTICLE_H

#include <iomanip>

#include "simprop/pid.h"
#include "simprop/units.h"
#include "simprop/utils/io.h"

namespace simprop {

class Particle {
  struct State {
    double z;
    double Gamma;
  };

  PID m_pid;
  State m_origin;
  State m_now;
  bool m_isPrimary;

 public:
  Particle(PID pid, double z, double Gamma, bool isPrimary = false)
      : m_pid(pid), m_origin({z, Gamma}), m_now({z, Gamma}), m_isPrimary(isPrimary) {}

  const State getNow() const { return m_now; }
  State& getNow() { return m_now; }
  const PID& getPid() const { return m_pid; }
  const State& getOrigin() const { return m_origin; }
  const double getRedshift() const { return m_now.z; }
  const double getGamma() const { return m_now.Gamma; }
  // const double getEnergy() const { return m_now.Gamma * getMassFromPid(m_pid); }
  const bool IsPrimary() const { return m_isPrimary; }

  friend std::ostream& operator<<(std::ostream& os, const Particle& p) {
    auto pidName = getPidName(p.m_pid);
    auto z = p.m_now.z;
    auto Gamma = p.m_now.Gamma;
    return os << pidName << " " << std::fixed << z << " " << std::scientific << Gamma;
  }
};

}  // namespace simprop

#endif