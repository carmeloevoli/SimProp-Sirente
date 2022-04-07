#include "simprop.h"

using namespace simprop;

int main() {
  try {
    utils::startup_information();
    auto particles = ParticleStack(proton, 1e6, 66);
    particles.buildInitialStates({0., 2.}, {1e8, 1e12}, 1.0);
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