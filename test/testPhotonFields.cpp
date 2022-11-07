
#include "gtest/gtest.h"
#include "simprop.h"

namespace simprop {

TEST(PhotonFields, CMB) {
  auto cmb = photonfields::CMB();
  EXPECT_DOUBLE_EQ(cmb.getMinPhotonEnergy(), 1e-5 * SI::eV);
  EXPECT_DOUBLE_EQ(cmb.getMaxPhotonEnergy(), 0.1 * SI::eV);
}

TEST(PhotonFields, Dominguez2011) {
  auto cmb = photonfields::CMB();
  EXPECT_DOUBLE_EQ(cmb.getMinPhotonEnergy(), 0.);
  EXPECT_DOUBLE_EQ(cmb.getMaxPhotonEnergy(), 0.);
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

}  // namespace simprop