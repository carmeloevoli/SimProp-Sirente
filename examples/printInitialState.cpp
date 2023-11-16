#include "simprop.h"

using namespace simprop;

void plot_initial_redshift() {
  auto cosmology = std::make_shared<cosmo::Planck2018>();
  RandomNumberGenerator rng = utils::RNG<double>(69);
  auto builder = SourceEvolutionBuilder(proton, {{1e10, 1e10}, {0., 3.0}, 2.2, 2.}, cosmology, 1e6);
  auto particles = builder.build(rng);
  utils::OutputFile out("SimProp_stack_redshift.txt");
  out << "#\n";
  for (const auto& p : particles) {
    out << p << "\n";
  }
}

void plot_initial_energy() {
  auto cosmology = std::make_shared<cosmo::Planck2018>();
  RandomNumberGenerator rng = utils::RNG<double>(96);
  auto builder = SourceEvolutionBuilder(proton, {{1e8, 1e13}, {0., 0.}, 2.2, 2.}, cosmology, 1e6);
  auto particles = builder.build(rng);
  utils::OutputFile out("SimProp_stack_energy.txt");
  out << "#\n";
  for (const auto& p : particles) {
    out << p << "\n";
  }
}

int main() {
  try {
    utils::startup_information();
    utils::Timer timer("main timer");
    plot_initial_redshift();
    plot_initial_energy();
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}