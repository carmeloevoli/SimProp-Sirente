#ifndef SIMPROP_PARTICLESTACK_H
#define SIMPROP_PARTICLESTACK_H

// #include <memory>
#include <vector>

#include "simprop/particle.h"
// #include "simprop/utils/random.h"

namespace simprop {

using ParticleStack = std::vector<Particle>;
using Range = std::pair<double, double>;

class Builder {
 protected:
  PID m_pid;
  size_t m_size;

 public:
  Builder(PID pid, size_t size = 1) : m_pid(pid), m_size(size) {}
  virtual ~Builder() = default;
  virtual ParticleStack build() const = 0;
};

class SingleParticleBuilder final : public Builder {
 protected:
  double m_z = 1;
  double m_Gamma = 1e12;

 public:
  SingleParticleBuilder(PID pid, size_t size = 1) : Builder(pid, size) {}
  void setRedshift(double z) { m_z = z; };
  void setGamma(double Gamma) { m_Gamma = Gamma; };

  ParticleStack build() const override {
    ParticleStack stack;
    stack.reserve(m_size);
    for (size_t i = 0; i < m_size; ++i) stack.emplace_back(Particle{m_pid, m_z, m_Gamma});
    assert(stack.size() == m_size);
    return stack;
  };
};

class SingleSourceBuilder final : public Builder {
 protected:
  double m_z = 1;
  Range m_GammaRange = {1e8, 1e14};
  double m_slope = 2.;

 public:
  SingleSourceBuilder(PID pid, size_t size = 1) : Builder(pid, size) {}
  void setRedshift(double z) { m_z = z; };
  void setGammaRange(Range GammaRange) { m_GammaRange = GammaRange; };
  void setSlope(double slope) { m_slope = slope; }

  ParticleStack build() const override {
    ParticleStack stack;
    // for (size_t i = 0; i < N; ++i) {
    //     auto z_i = getRndRedshift(zRange, 2, m_rng());
    //     auto Gamma_i = getRndGamma(gammaRange, slope, m_rng());
    //     m_particles.push_back(Particle{m_pid, z_i, Gamma_i});
    //   }
    //   LOGD << "built primaries with size " << m_particles.size();
    //   auto z_r = getRedshiftRange();
    //   LOGD << "z range (" << z_r.first << "," << z_r.second << ")";
    //   auto G_r = getGammaRange();
    //   LOGD << "Gamma range (" << G_r.first << "," << G_r.second << ")";
    //}}
  };

  // ParticleStack buildSingleEnergyParticleStack(PID pid, double Gamma, double z, size_t N = 1);
  // ParticleStack buildSingleSourceParticleStack(PID pid, Range GammaRange, double z, double slope,
  //                                              size_t N);

  // class ParticleStack {
  //  public:
  //   using stack = std::vector<Particle>;
  //   using iterator = typename stack::iterator;
  //   using const_iterator = typename stack::const_iterator;

  //  protected:
  //   PID m_pid;
  //   size_t m_size;
  //   stack m_particles;
  //   RandomNumberGenerator& m_rng;

  //  public:
  //   iterator begin() { return m_particles.begin(); }
  //   iterator end() { return m_particles.end(); }
  //   const_iterator begin() const { return m_particles.begin(); }
  //   const_iterator end() const { return m_particles.end(); }

  //  public:
  //   explicit ParticleStack(PID pid, int nParticles, RandomNumberGenerator& rng);
  //   virtual ~ParticleStack() = default;

  //   void buildInitialState(Range zRange, Range gammaRange, double slope);
  //   void buildSingleParticleStack(double z, double Gamma);
  //   void buildMultipleParticleStack(Range zRange, Range gammaRange, double slope);

  //   Range getRedshiftRange() const;
  //   Range getGammaRange() const;
  // };

}  // namespace simprop

#endif  // SIMPROP_PARTICLESTACK_H