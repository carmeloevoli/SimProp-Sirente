void decay() {
  //   int charge = input->GetCharge();
  //   double E = input->GetEprod();
  //   double z = input->GetZprod();
  //   int nBr = input->GetBranch() + 1;
  //   assert(E >= 0. && z >= 0.);
  //   double gamma = E / mmu;
  //   double E1, E2, E3, p1, p2, p3;  // quantities in CoM frame
  //   // 1 = muon neutrino, 2 = electron neutrino, 3 = electron
  //   const double Enumax = (mmu - me * me / mmu) / 2.;
  //   do {
  //     E1 = gRandom->Uniform(Enumax);
  //     E2 = gRandom->Uniform(Enumax);
  //     if (E1 + E2 < Enumax) {
  //       E1 = Enumax - E1;
  //       E2 = Enumax - E2;
  //     }  // (E1, E2) uniform in triangle (0, Enumax), (Enumax, Enumax), (Enumax, 0)
  //     E3 = mmu - E1 - E2;
  //     p1 = E1;
  //     p2 = E2;
  //     p3 = sqrt(E3 * E3 - me * me);
  //   } while (E3 < me || p3 > p1 + p2 || p1 > p2 + p3 || p2 > p3 + p1);

  //   double c12 = (p3 * p3 - p1 * p1 - p2 * p2) / (2. * p1 * p2);  // angle btw neutrinos
  //   double s12 = sqrt(1. - c12 * c12);
  //   double c1 = gRandom->Uniform(-1., 1.);  // angle btw neutrino 1 and line of sight
  //   double s1 = sqrt(1. - c1 * c1);
  //   double ca = cos(gRandom->Uniform(6.283185307));  // angle bw nu2 and (nu1, line of sight)
  //   plane double c2 = c12 * c1 - s12 * s1 * ca;

  //   double Enu1 = gamma * (E1 + p1 * c1);  // in lab frame
  //   double Enu2 = gamma * (E2 + p2 * c2);
  //   double Ee = E - Enu1 - Enu2;

  //   input->SetEint(E);
  //   input->SetZint(z);
  //   input->SetIntMult(1003);
  //   vector<Particle> output;
  //   output.push_back(Particle(0, 0, nBr, Enu1, z, (charge > 0) ? eAntineutrino_mu :
  //   eNeutrino_mu)); output.push_back(Particle(0, 0, nBr, Enu2, z, (charge > 0) ? eNeutrino_e :
  //   eAntineutrino_e)); output.push_back(Particle(0, charge, nBr, Ee, z, eElectron)); cout <<
  //   "electron " << Ee << "; neutrino " << Enu1 << "; neutrino " << Enu2 << endl;
}