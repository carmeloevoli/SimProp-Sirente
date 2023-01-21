std::shared_ptr<interactions::PhotoPionProduction> m_pppcmb;

void buildStochasticInteractions() {
  auto sigma = std::make_shared<xsecs::PhotoPionProductionXsec>();
  m_pppcmb = std::make_shared<interactions::PhotoPionProduction>(sigma, m_cmb);
}

double computeStochasticRedshiftInterval(const Particle& particle) {
  const auto pid = particle.getPid();
  const auto Gamma = particle.getGamma();
  const auto zNow = particle.getRedshift();
  const auto dtdz = m_cosmology->dtdz(zNow);
  const auto dt = std::fabs(1. / m_pppcmb->rate(pid, Gamma, zNow));
  // TODO why to put the fabs?
  return -dt / dtdz * std::log(1. - m_rng());
}

void run(std::string filename) {
  assert(m_pppcmb != nullptr);
  auto it = m_stack.begin();
  size_t iniSize = m_stack.size();
  size_t counter = 0;
  const double initEnergy = sumEnergy();
  while (it != m_stack.end()) {
    if (counter % iniSize == 0) {
      LOGD << counter / iniSize << "\t" << sumEnergy() / initEnergy << "\t" << m_stack.size();
    }

    const auto distance = it - m_stack.begin();
    const auto nowRedshift = it->getRedshift();
    const auto dz_s = computeStochasticRedshiftInterval(*it);
    const auto dz_c = computeLossesRedshiftInterval(*it);
    assert(dz_s > 0. && dz_c > 0. && dz_c <= nowRedshift);

    if (dz_s > dz_c || dz_s > nowRedshift) {
      const auto Gamma = it->getGamma();
      const auto dz = dz_c;
      const auto deltaGamma = computeDeltaGamma(*it, dz);
      it->getNow() = {nowRedshift - dz, Gamma * (1. - deltaGamma)};
    } else {
      it->deactivate();
      const auto dz = dz_s;
      auto finalState = m_pppcmb->finalState(*it, nowRedshift - dz, m_rng);
      m_stack.insert(m_stack.end(), finalState.begin(), finalState.end());
    }
    it = std::find_if(m_stack.begin() + distance, m_stack.end(), IsActive);
    counter++;
  }
}  // run()

void evolvePopulation() {
  for (size_t i = 100; i < 1000; ++i) {
    RandomNumberGenerator rng = utils::RNG<double>(Seed(i));
    utils::Timer timer("timer for first test");
    Evolutor evolutor(rng);
    evolutor.buildCosmologicalParticleStack(2.7, 0., Redshift(1.0), 1e5);
    evolutor.buildPhotonFields();
    evolutor.buildContinuousLosses();
    evolutor.run("test_proton_cosmology.txt");
    evolutor.dump("test_" + std::to_string(i) + "_spectrum_z1.0_m0.txt");
  }
}
