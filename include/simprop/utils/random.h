#ifndef SIMPROP_UTILS_RANDOM_H
#define SIMPROP_UTILS_RANDOM_H

#include <random>

namespace simprop {
namespace utils {

template <class FloatType = double,
          class = std::enable_if_t<std::is_floating_point<FloatType>::value> >
class RNG {
 public:
  typedef FloatType result_type;
  typedef std::mt19937_64 generator_type;
  typedef std::uniform_real_distribution<FloatType> distribution_type;

  explicit RNG(const unsigned int seed) { eng = generator_type(seed); }

  // generate next random value in distribution
  result_type operator()() { return dist(eng); }
  // will always yield 0.0 for this class type
  constexpr result_type min() const { return dist.min(); }
  // will always yield 1.0 for this class type
  constexpr result_type max() const { return dist.max(); }
  // resets internal state such that next call to operator()
  // does not rely on previous call
  void reset_distribution_state() { dist.reset(); }
  // uniform distribution
  result_type uniform(double vMin, double vMax) { return dist(eng) * (vMax - vMin) + vMin; }

 private:
  generator_type eng;
  distribution_type dist;
};

}  // namespace utils

using RandomNumberGenerator = simprop::utils::RNG<double>;

}  // namespace simprop

#endif  // SIMPROP_UTILS_RANDOM_H