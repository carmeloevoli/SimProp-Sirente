#ifndef SIMPROP_INTERACTIONS_PHOTODISINTEGRATION_H
#define SIMPROP_INTERACTIONS_PHOTODISINTEGRATION_H

namespace simprop {
namespace interactions {

class PhotoDisintegration final : public Interaction {
 public:
  PhotoDisintegration(const std::shared_ptr<photonfields::PhotonField>& phField);
  virtual ~PhotoDisintegration() = default;
  double rate(PID pid, double Gamma, double z = 0) const override;

  std::vector<Particle> finalState(const Particle& particle, double zInteractionPoint,
                                   RandomNumberGenerator& rng) const override;
};

}  // namespace simprop

#endif  // SIMPROP_INTERACTIONS_PHOTODISINTEGRATION_H