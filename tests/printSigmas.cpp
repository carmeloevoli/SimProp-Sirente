#include "simprop.h"

using namespace simprop;

int main() {
  auto photoPion = interactions::PhotoPion();
  auto energyAxis = utils::LogAxis(0.1 * SI::GeV, 1e12 * SI::GeV, 1000);
  for (auto E : energyAxis) {
    std::cout << E / SI::GeV << "\t";
    std::cout << photoPion.sigma(E) / SI::mbarn << "\n";
  }
}