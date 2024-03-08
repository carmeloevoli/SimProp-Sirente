#include "simprop/interactions/PhotoPionProductionSophia.h"

#include "simprop/utils/logging.h"
#include "sophia_interface.h"

namespace simprop {
namespace interactions {

PhotoPionProductionSophia::PhotoPionProductionSophia(
    const std::shared_ptr<photonfields::PhotonField>& phField)
    : PhotoPionProduction(phField) {
  LOGD << "calling " << __func__ << " constructor";
}

/*
list of output particle IDs
   13  proton
  -13  anti-proton
   14  neutron
  -14  anti-neutron
    1   photon
    2   e+
    3   e-
   15  nu_e
   16  anti-nu_e
   17  nu_mu
   18  anti-nu_mu
*/
PID ID_sophia_to_SimProp(int sophiaID) {
  switch (sophiaID) {
    case 13:
      return proton;
    case 14:
      return proton;  // TODO neutron -> proton?
    case -13:
      return antiproton;
    case -14:  // TODO antineutron == neutron?
      return neutron;
    case 1:
      return photon;
    case 2:
      return positron;
    case 3:
      return electron;
    case 15:
      return neutrino_e;
    case 16:
      return antineutrino_e;
    case 17:
      return neutrino_mu;
    case 18:
      return antineutrino_mu;
    case 6:
      return pionNeutral;
    case 7:
      return pionPlus;
    case 8:
      return pionMinus;
    default:
      throw std::runtime_error("ID_sophia_to_SimProp: unkown particle ID");
  }
}

std::vector<Particle> PhotoPionProductionSophia::finalState(const Particle& incomingParticle,
                                                            double zInteractionPoint,
                                                            RandomNumberGenerator& rng) const {
  const auto pid = incomingParticle.getPid();
  assert(pidIsNucleus(pid));
  const auto w = incomingParticle.getWeight();
  const auto Gamma = incomingParticle.getGamma();

  const auto nucleon = pidIsNucleon(pid) ? pid : pickNucleon(rng(), pid);
  const auto nucleonEnergy = Gamma * SI::protonMassC2;
  const auto photonEnergy = sampleEps(rng(), nucleon, nucleonEnergy, zInteractionPoint);

  const bool onProton = (nucleon == proton);
  const double Ein = nucleonEnergy / SI::GeV;
  const double eps = photonEnergy / SI::GeV;
  const bool declareChargedPionsStable = true;

  sophia_interface SI;
  sophiaevent_output seo = SI.sophiaevent(onProton, Ein, eps, declareChargedPionsStable);

  std::vector<Particle> outgoingParticle;
  int Nout = seo.Nout;
  for (int i = 0; i < Nout; ++i) {
    PID pid = ID_sophia_to_SimProp(seo.outPartID[i]);
    auto mass = getPidMass(pid);
    auto E = seo.outPartP[3][i] * SI::GeV;
    if (mass > 0.) {
      outgoingParticle.push_back({pid, zInteractionPoint, E / mass, w});
    } else {
      outgoingParticle.push_back({pid, zInteractionPoint, E, w});  // TODO terrible patch :(
    }
  }
  return outgoingParticle;
}

}  // namespace interactions
}  // namespace simprop