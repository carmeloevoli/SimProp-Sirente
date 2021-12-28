#include "simprop.h"

using namespace simprop;

int main() {
  try {
    utils::startup_information();
    auto particles = ParticleStack(proton, 10000, 333);
    particles.buildInitialStates({0., 2.}, {1e17 * SI::eV, 1e21 * SI::eV}, 1.0);
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