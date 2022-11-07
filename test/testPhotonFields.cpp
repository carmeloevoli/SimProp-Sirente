
#include "gtest/gtest.h"
#include "simprop.h"

namespace simprop {

TEST(PhotonFields, CMB) {
  auto cmb = photonfields::CMB();
  EXPECT_DOUBLE_EQ(cmb.getMinPhotonEnergy(), 1e-5 * SI::eV);
  EXPECT_DOUBLE_EQ(cmb.getMaxPhotonEnergy(), 0.1 * SI::eV);
  EXPECT_DOUBLE_EQ(cmb.density(SI::eV), 0.);
  EXPECT_DOUBLE_EQ(cmb.density(1e-6 * SI::eV), 0.);
  EXPECT_DOUBLE_EQ(cmb.I_gamma(SI::eV), 0.);
  EXPECT_DOUBLE_EQ(cmb.I_gamma(1e-6 * SI::eV), 0.);
}

TEST(PhotonFields, Dominguez2011) {
  auto ebl = photonfields::Dominguez2011PhotonField();
  EXPECT_NEAR(ebl.getMinPhotonEnergy(), 0.0012388 * SI::eV, 1e-6 * SI::eV);
  EXPECT_NEAR(ebl.getMaxPhotonEnergy(), 12.2744 * SI::eV, 1e-3 * SI::eV);
  EXPECT_DOUBLE_EQ(ebl.density(1e-3 * SI::eV), 0.);
  EXPECT_DOUBLE_EQ(ebl.density(1e2 * SI::eV), 0.);
  EXPECT_DOUBLE_EQ(ebl.I_gamma(1e-3 * SI::eV), 0.);
  EXPECT_DOUBLE_EQ(ebl.I_gamma(1e2 * SI::eV), 0.);
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

}  // namespace simprop