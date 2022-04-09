#include <memory>

#include "gtest/gtest.h"
#include "simprop.h"

namespace simprop {

TEST(Pid, elementaryName) {
  EXPECT_STREQ("photon", getPidName(photon).c_str());
  EXPECT_STREQ("nu_e", getPidName(neutrino_e).c_str());
  EXPECT_STREQ("antinu_e", getPidName(antineutrino_e).c_str());
  EXPECT_STREQ("nu_mu", getPidName(neutrino_mu).c_str());
  EXPECT_STREQ("antinu_mu", getPidName(antineutrino_mu).c_str());
  EXPECT_STREQ("electron", getPidName(electron).c_str());
  EXPECT_STREQ("positron", getPidName(positron).c_str());
  EXPECT_STREQ("pion_0", getPidName(pionNeutral).c_str());
  EXPECT_STREQ("pion_plus", getPidName(pionPlus).c_str());
  EXPECT_STREQ("pion_minus", getPidName(pionMinus).c_str());
}

TEST(Pid, elementaryNotNucleus) {
  EXPECT_FALSE(isNucleus(photon));
  EXPECT_FALSE(isNucleus(neutrino_e));
  EXPECT_FALSE(isNucleus(antineutrino_e));
  EXPECT_FALSE(isNucleus(neutrino_mu));
  EXPECT_FALSE(isNucleus(antineutrino_mu));
  EXPECT_FALSE(isNucleus(electron));
  EXPECT_FALSE(isNucleus(positron));
  EXPECT_FALSE(isNucleus(pionNeutral));
  EXPECT_FALSE(isNucleus(pionPlus));
  EXPECT_FALSE(isNucleus(pionMinus));
  EXPECT_TRUE(isNucleus(proton));
  EXPECT_TRUE(isNucleus(neutron));
  EXPECT_TRUE(isNucleus(antiproton));
  EXPECT_TRUE(isNucleus(deuterium));
}

TEST(Pid, nucleiName) {
  EXPECT_STREQ("proton", getPidName(proton).c_str());
  EXPECT_STREQ("neutron", getPidName(neutron).c_str());
  EXPECT_STREQ("antiproton", getPidName(antiproton).c_str());
  EXPECT_STREQ("deuterium", getPidName(deuterium).c_str());
  EXPECT_STREQ("He3", getPidName(He3).c_str());
  EXPECT_STREQ("He4", getPidName(He4).c_str());
  EXPECT_STREQ("Be9", getPidName(Be9).c_str());
  EXPECT_STREQ("B10", getPidName(B10).c_str());
  EXPECT_STREQ("B11", getPidName(B11).c_str());
  EXPECT_STREQ("C12", getPidName(C12).c_str());
  EXPECT_STREQ("C13", getPidName(C13).c_str());
  EXPECT_STREQ("N14", getPidName(N14).c_str());
  EXPECT_STREQ("N15", getPidName(N15).c_str());
  EXPECT_STREQ("O16", getPidName(O16).c_str());
  EXPECT_STREQ("O17", getPidName(O17).c_str());
  EXPECT_STREQ("O18", getPidName(O18).c_str());
  EXPECT_STREQ("F19", getPidName(F19).c_str());
  EXPECT_STREQ("Ne20", getPidName(Ne20).c_str());
  EXPECT_STREQ("Ne21", getPidName(Ne21).c_str());
  EXPECT_STREQ("Ne22", getPidName(Ne22).c_str());
  EXPECT_STREQ("Na23", getPidName(Na23).c_str());
  EXPECT_STREQ("Mg24", getPidName(Mg24).c_str());
  EXPECT_STREQ("Mg25", getPidName(Mg25).c_str());
  EXPECT_STREQ("Mg26", getPidName(Mg26).c_str());
  EXPECT_STREQ("Al27", getPidName(Al27).c_str());
  EXPECT_STREQ("Si28", getPidName(Si28).c_str());
  EXPECT_STREQ("Si29", getPidName(Si29).c_str());
  EXPECT_STREQ("Si30", getPidName(Si30).c_str());
  EXPECT_STREQ("P31", getPidName(P31).c_str());
  EXPECT_STREQ("S32", getPidName(S32).c_str());
  EXPECT_STREQ("S33", getPidName(S33).c_str());
  EXPECT_STREQ("S34", getPidName(S34).c_str());
  EXPECT_STREQ("Cl35", getPidName(Cl35).c_str());
  EXPECT_STREQ("Ar36", getPidName(Ar36).c_str());
  EXPECT_STREQ("Cl37", getPidName(Cl37).c_str());
  EXPECT_STREQ("Ar38", getPidName(Ar38).c_str());
  EXPECT_STREQ("K39", getPidName(K39).c_str());
  EXPECT_STREQ("Ca40", getPidName(Ca40).c_str());
  EXPECT_STREQ("K41", getPidName(K41).c_str());
  EXPECT_STREQ("Ca42", getPidName(Ca42).c_str());
  EXPECT_STREQ("Ca43", getPidName(Ca43).c_str());
  EXPECT_STREQ("Ca44", getPidName(Ca44).c_str());
  EXPECT_STREQ("Sc45", getPidName(Sc45).c_str());
  EXPECT_STREQ("Ti46", getPidName(Ti46).c_str());
  EXPECT_STREQ("Ti47", getPidName(Ti47).c_str());
  EXPECT_STREQ("Ti48", getPidName(Ti48).c_str());
  EXPECT_STREQ("Ti49", getPidName(Ti49).c_str());
  EXPECT_STREQ("Cr50", getPidName(Cr50).c_str());
  EXPECT_STREQ("V51", getPidName(V51).c_str());
  EXPECT_STREQ("Cr52", getPidName(Cr52).c_str());
  EXPECT_STREQ("Cr53", getPidName(Cr53).c_str());
  EXPECT_STREQ("Fe54", getPidName(Fe54).c_str());
  EXPECT_STREQ("Mn55", getPidName(Mn55).c_str());
  EXPECT_STREQ("Fe56", getPidName(Fe56).c_str());
}

TEST(Pid, nucleiMass) {
  EXPECT_EQ(1, getNucleusMassNumber(proton));
  EXPECT_EQ(1, getNucleusMassNumber(neutron));
  EXPECT_EQ(1, getNucleusMassNumber(antiproton));
  EXPECT_EQ(2, getNucleusMassNumber(deuterium));
  EXPECT_EQ(3, getNucleusMassNumber(He3));
  EXPECT_EQ(4, getNucleusMassNumber(He4));
  EXPECT_EQ(9, getNucleusMassNumber(Be9));
  EXPECT_EQ(10, getNucleusMassNumber(B10));
  EXPECT_EQ(11, getNucleusMassNumber(B11));
  EXPECT_EQ(12, getNucleusMassNumber(C12));
  EXPECT_EQ(13, getNucleusMassNumber(C13));
  EXPECT_EQ(14, getNucleusMassNumber(N14));
  EXPECT_EQ(15, getNucleusMassNumber(N15));
  EXPECT_EQ(16, getNucleusMassNumber(O16));
  EXPECT_EQ(17, getNucleusMassNumber(O17));
  EXPECT_EQ(18, getNucleusMassNumber(O18));
  EXPECT_EQ(19, getNucleusMassNumber(F19));
  EXPECT_EQ(20, getNucleusMassNumber(Ne20));
  EXPECT_EQ(21, getNucleusMassNumber(Ne21));
  EXPECT_EQ(22, getNucleusMassNumber(Ne22));
  EXPECT_EQ(23, getNucleusMassNumber(Na23));
  EXPECT_EQ(24, getNucleusMassNumber(Mg24));
  EXPECT_EQ(25, getNucleusMassNumber(Mg25));
  EXPECT_EQ(26, getNucleusMassNumber(Mg26));
  EXPECT_EQ(27, getNucleusMassNumber(Al27));
  EXPECT_EQ(28, getNucleusMassNumber(Si28));
  EXPECT_EQ(29, getNucleusMassNumber(Si29));
  EXPECT_EQ(30, getNucleusMassNumber(Si30));
  EXPECT_EQ(31, getNucleusMassNumber(P31));
  EXPECT_EQ(32, getNucleusMassNumber(S32));
  EXPECT_EQ(33, getNucleusMassNumber(S33));
  EXPECT_EQ(34, getNucleusMassNumber(S34));
  EXPECT_EQ(35, getNucleusMassNumber(Cl35));
  EXPECT_EQ(36, getNucleusMassNumber(Ar36));
  EXPECT_EQ(37, getNucleusMassNumber(Cl37));
  EXPECT_EQ(38, getNucleusMassNumber(Ar38));
  EXPECT_EQ(39, getNucleusMassNumber(K39));
  EXPECT_EQ(40, getNucleusMassNumber(Ca40));
  EXPECT_EQ(41, getNucleusMassNumber(K41));
  EXPECT_EQ(42, getNucleusMassNumber(Ca42));
  EXPECT_EQ(43, getNucleusMassNumber(Ca43));
  EXPECT_EQ(44, getNucleusMassNumber(Ca44));
  EXPECT_EQ(45, getNucleusMassNumber(Sc45));
  EXPECT_EQ(46, getNucleusMassNumber(Ti46));
  EXPECT_EQ(47, getNucleusMassNumber(Ti47));
  EXPECT_EQ(48, getNucleusMassNumber(Ti48));
  EXPECT_EQ(49, getNucleusMassNumber(Ti49));
  EXPECT_EQ(50, getNucleusMassNumber(Cr50));
  EXPECT_EQ(51, getNucleusMassNumber(V51));
  EXPECT_EQ(52, getNucleusMassNumber(Cr52));
  EXPECT_EQ(53, getNucleusMassNumber(Cr53));
  EXPECT_EQ(54, getNucleusMassNumber(Fe54));
  EXPECT_EQ(55, getNucleusMassNumber(Mn55));
  EXPECT_EQ(56, getNucleusMassNumber(Fe56));
}

TEST(Pid, nucleiCharge) {
  EXPECT_EQ(1, getNucleusCharge(proton));
  EXPECT_EQ(0, getNucleusCharge(neutron));
  EXPECT_EQ(-1, getNucleusCharge(antiproton));
  EXPECT_EQ(1, getNucleusCharge(deuterium));
  EXPECT_EQ(2, getNucleusCharge(He3));
  EXPECT_EQ(2, getNucleusCharge(He4));
  EXPECT_EQ(4, getNucleusCharge(Be9));
  EXPECT_EQ(5, getNucleusCharge(B10));
  EXPECT_EQ(5, getNucleusCharge(B11));
  EXPECT_EQ(6, getNucleusCharge(C12));
  EXPECT_EQ(6, getNucleusCharge(C13));
  EXPECT_EQ(7, getNucleusCharge(N14));
  EXPECT_EQ(7, getNucleusCharge(N15));
  EXPECT_EQ(8, getNucleusCharge(O16));
  EXPECT_EQ(8, getNucleusCharge(O17));
  EXPECT_EQ(8, getNucleusCharge(O18));
  EXPECT_EQ(9, getNucleusCharge(F19));
  EXPECT_EQ(10, getNucleusCharge(Ne20));
  EXPECT_EQ(10, getNucleusCharge(Ne21));
  EXPECT_EQ(10, getNucleusCharge(Ne22));
  EXPECT_EQ(11, getNucleusCharge(Na23));
  EXPECT_EQ(12, getNucleusCharge(Mg24));
  EXPECT_EQ(12, getNucleusCharge(Mg25));
  EXPECT_EQ(12, getNucleusCharge(Mg26));
  EXPECT_EQ(13, getNucleusCharge(Al27));
  EXPECT_EQ(14, getNucleusCharge(Si28));
  EXPECT_EQ(14, getNucleusCharge(Si29));
  EXPECT_EQ(14, getNucleusCharge(Si30));
  EXPECT_EQ(15, getNucleusCharge(P31));
  EXPECT_EQ(16, getNucleusCharge(S32));
  EXPECT_EQ(16, getNucleusCharge(S33));
  EXPECT_EQ(16, getNucleusCharge(S34));
  EXPECT_EQ(17, getNucleusCharge(Cl35));
  EXPECT_EQ(18, getNucleusCharge(Ar36));
  EXPECT_EQ(17, getNucleusCharge(Cl37));
  EXPECT_EQ(18, getNucleusCharge(Ar38));
  EXPECT_EQ(19, getNucleusCharge(K39));
  EXPECT_EQ(20, getNucleusCharge(Ca40));
  EXPECT_EQ(19, getNucleusCharge(K41));
  EXPECT_EQ(20, getNucleusCharge(Ca42));
  EXPECT_EQ(20, getNucleusCharge(Ca43));
  EXPECT_EQ(20, getNucleusCharge(Ca44));
  EXPECT_EQ(21, getNucleusCharge(Sc45));
  EXPECT_EQ(22, getNucleusCharge(Ti46));
  EXPECT_EQ(22, getNucleusCharge(Ti47));
  EXPECT_EQ(22, getNucleusCharge(Ti48));
  EXPECT_EQ(22, getNucleusCharge(Ti49));
  EXPECT_EQ(24, getNucleusCharge(Cr50));
  EXPECT_EQ(23, getNucleusCharge(V51));
  EXPECT_EQ(24, getNucleusCharge(Cr52));
  EXPECT_EQ(24, getNucleusCharge(Cr53));
  EXPECT_EQ(26, getNucleusCharge(Fe54));
  EXPECT_EQ(25, getNucleusCharge(Mn55));
  EXPECT_EQ(26, getNucleusCharge(Fe56));
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

}  // namespace simprop