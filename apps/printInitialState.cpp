#include "simprop.h"

using namespace simprop;

int main() {
  try {
    utils::startup_information();
    utils::Timer timer("main timer");
    RandomNumberGenerator rng = utils::RNG<double>(Seed(69));
    {
      auto builder = SingleSourceBuilder(proton, {{1e8, 1e12}, 1.0, 1.0, 1e12}, 1e6);
      auto particles = builder.build(rng);
      utils::OutputFile out("test_singlesource_initial.txt");
      out << "#\n";
      for (const auto& p : particles) {
        out << p << "\n";
      }
    }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}