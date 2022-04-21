#include "simprop.h"

using namespace simprop;

int main() {
  try {
    utils::startup_information();
    utils::Timer timer("main timer");
    RandomNumberGenerator rng = utils::RNG<double>(66);
    auto builder = SourceEvolutionBuilder(proton, 1e6);
    builder.setGammaRange({1e8, 1e12});
    builder.setRedshiftRange({0., 2.});
    builder.setSlope(1.0);
    builder.setEvolutionIndex(3.0);
    auto particles = builder.build(rng);
    utils::OutputFile out("test_initial.txt");
    out << "#\n";
    for (auto& p : particles) {
      out << p << "\n";
    }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}