#include "simprop.h"

using namespace simprop;

void plot_photonfields_local(const photonfields::PhotonFields& phFields) {
  const auto ePhoton = utils::LogAxis<double>(1e-5 * SI::eV, 1e2 * SI::eV, 1000);
  utils::OutputFile out("SimProp_photonfields_z0.txt");
  const auto units = SI::nW / pow2(SI::meter) / SI::sr / SI::cOver4pi;
  out << "# energy [eV] - lambda [mu] - density []\n";
  out << std::scientific;
  for (auto E : ePhoton) {
    out << E / SI::eV << "\t";
    out << energyToWavelenght(E) / SI::micron << "\t";
    for (auto phField : phFields) out << pow2(E) * phField->density(E) / units << "\t";
    out << "\n";
  }
}

void plot_photonfields_Igamma(const photonfields::PhotonFields& phFields) {
  const auto ePhoton = utils::LogAxis<double>(1e-5 * SI::eV, 1e2 * SI::eV, 1000);
  utils::OutputFile out("SimProp_photonfields_Igamma.txt");
  const auto units = 1. / SI::eV / SI::m3;
  out << "# energy [eV] - wavelength [mu] - I_gamma [eV^-1 m^-3] - \n";
  out << std::scientific;
  for (auto E : ePhoton) {
    out << E / SI::eV << "\t";
    out << energyToWavelenght(E) / SI::micron << "\t";
    for (auto phField : phFields) out << E * phField->I_gamma(E) / units << "\t";
    out << "\n";
  }
}

void plot_photonfields_density(const photonfields::PhotonFields& phFields) {
  const auto redshifts = utils::LinAxis<double>(0, 6, 1000);
  utils::OutputFile out("SimProp_photonfields_density.txt");
  out << "# z - photon density\n";
  out << std::scientific;
  auto units = 1. / SI::m3;
  for (auto z : redshifts) {
    out << z << "\t";
    for (auto phField : phFields)
      out << phField->computeIntegratedDensity(z) / std::pow(1. + z, 3.) / units << "\t";
    out << "\n";
  }
}

void plot_photonfields_redshifts(const photonfields::PhotonFields& phFields) {
  const auto units = SI::nW / pow2(SI::meter) / SI::sr / SI::cOver4pi;
  const auto cmb = std::make_shared<photonfields::CMB>();
  utils::OutputFile out("SimProp_photonfields_redshifts.txt");
  out << "# eps [eV] - photon density\n";
  out << std::scientific;
  const auto ePhoton = utils::LogAxis<double>(1e-5 * SI::eV, 1e2 * SI::eV, 1000);
  std::vector<double> redshifts{0., 0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0};
  for (auto E : ePhoton) {
    out << E / SI::eV << "\t";
    out << energyToWavelenght(E) / SI::micron << "\t";
    for (auto z : redshifts) {
      for (auto phField : phFields)
        out << pow2(E) * phField->density(E, z) / std::pow(1. + z, 3.) / units << "\t";
    }
    out << "\n";
  }
}

void plot_photonfields_radio() {
  const auto units = SI::watt / SI::Hz / pow2(SI::meter) / SI::sr;
  const auto radio = std::make_shared<photonfields::Nitu2021RadioPhotonField>();
  utils::OutputFile out("SimProp_photonfields_radio.txt");
  out << "# nu [Hz] - photon density\n";
  out << std::scientific;
  const auto nuAxis = utils::LogAxis<double>(1e3 * SI::Hz, 1e11 * SI::Hz, 1000);
  std::vector<double> redshifts{0., 0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0};
  for (auto nu : nuAxis) {
    out << nu / SI::Hz << "\t";
    auto E = nu * SI::hPlanck;
    auto factor = SI::cLight / 4. * M_PI * (SI::hPlanck * E);
    out << E / SI::eV << "\t";
    for (auto z : redshifts) {
      out << factor * radio->density(E, z) / std::pow(1. + z, 3.) / units << "\t";
    }
    out << "\n";
  }
}

int main() {
  try {
    utils::startup_information();
    // std::vector<std::shared_ptr<photonfields::PhotonField> > phFields;
    photonfields::PhotonFields phFields;

    // phFields.push_back(std::make_shared<photonfields::CMB>());
    // phFields.push_back(std::make_shared<photonfields::Dominguez2011PhotonField>());
    // phFields.push_back(std::make_shared<photonfields::Dominguez2011LowerPhotonField>());
    // phFields.push_back(std::make_shared<photonfields::Dominguez2011UpperPhotonField>());
    // phFields.push_back(std::make_shared<photonfields::Saldana2021PhotonField>());
    // phFields.push_back(std::make_shared<photonfields::Saldana2021LowerPhotonField>());
    // phFields.push_back(std::make_shared<photonfields::Saldana2021UpperPhotonField>());
    // phFields.push_back(std::make_shared<photonfields::Gilmore2012PhotonField>());

    // plot_photonfields_local(phFields);
    // plot_photonfields_Igamma(phFields);
    // plot_photonfields_density(phFields);
    // plot_photonfields_redshifts(phFields);
    plot_photonfields_radio();
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}
