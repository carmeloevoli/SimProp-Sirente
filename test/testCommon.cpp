#include <memory>

#include "gtest/gtest.h"
#include "simprop.h"

namespace simprop {

TEST(Common, energyToFrequency) {
  double photonEnergy = SI::eV;
  double photonFrenquency = energyToFrequency(photonEnergy);
  EXPECT_NEAR(photonFrenquency, 2.418e14 * SI::Hz, 1e11 * SI::Hz);
}

TEST(Common, energyToWavelenght) {
  double photonEnergy = SI::eV;
  double photonWavelenght = energyToWavelenght(photonEnergy);
  EXPECT_NEAR(photonWavelenght, 1.23984193 * SI::micron, 1e-3 * SI::micron);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

}  // namespace simprop