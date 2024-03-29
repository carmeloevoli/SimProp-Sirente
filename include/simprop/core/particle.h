#ifndef SIMPROP_PARTICLE_H
#define SIMPROP_PARTICLE_H

#include <iomanip>

#include "simprop/core/pid.h"
#include "simprop/core/units.h"
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
  double m_weight;
  bool m_doPropagate;

 public:
  Particle(PID pid, double z, double Gamma, double weight = 1.)
      : m_pid(pid),
        m_origin({z, Gamma}),
        m_now({z, Gamma}),
        m_weight(weight),
        m_doPropagate(true) {}

  // Copy constructor
  Particle(const Particle& particle) {
    m_pid = particle.m_pid;
    m_origin = particle.m_origin;
    m_now = particle.m_now;
    m_weight = particle.m_weight;
    m_doPropagate = particle.m_doPropagate;
  }

  const State getNow() const { return m_now; }
  State& getNow() { return m_now; }
  const PID& getPid() const { return m_pid; }
  const State& getOrigin() const { return m_origin; }
  const double getRedshift() const { return m_now.z; }
  const double getGamma() const { return m_now.Gamma; }
  const double getEnergy() const { return m_now.Gamma * getPidMass(m_pid); }
  const double getWeight() const { return m_weight; }
  const bool isNucleus() const { return pidIsNucleus(m_pid); }
  const bool isActive() const { return m_doPropagate; }

  void deactivate() { m_doPropagate = false; }
  void activate() { m_doPropagate = true; }

  friend std::ostream& operator<<(std::ostream& os, const Particle& p) {
    auto pidName = getPidName(p.m_pid);
    auto z = p.m_now.z;
    auto Gamma = p.m_now.Gamma;
    auto weight = p.m_weight;
    auto doPropagate = p.m_doPropagate;
    return os << pidName << " " << std::scientific << z << " " << Gamma << " " << weight << " "
              << std::boolalpha << doPropagate;
  }
};

using ParticleStack = std::vector<Particle>;

}  // namespace simprop

#endif