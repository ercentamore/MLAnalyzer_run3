#include "MLAnalyzer_run3/RecHitAnalyzer/interface/RecHitAnalyzer.h"
#include <algorithm>

using std::vector;
using namespace trigger;

std::vector<unsigned int> vAs_ditau;
std::vector<unsigned int> vAs_diphoton;
std::vector<unsigned int> vPhotons;

std::vector<unsigned int> vReco_First_Photons_Idxs;
std::vector<unsigned int> vReco_Second_Photons_Idxs;
std::vector<unsigned int> vGen_As_diphoton_Idxs;
std::vector<unsigned int> vGen_As_ditau_Idxs;

std::vector<float> vA_diphoton_gen_m0_;
std::vector<float> vA_diphoton_gen_dR_;
std::vector<float> vA_diphoton_gen_E_;
std::vector<float> vA_diphoton_gen_pT_;
std::vector<float> vA_diphoton_gen_eta_;
std::vector<float> vA_diphoton_gen_phi_;
std::vector<float> vA_diphoton_reco_M_;
std::vector<float> vA_diphoton_reco_dR_;
std::vector<float> vA_diphoton_reco_E_;
std::vector<float> vA_diphoton_reco_pT_;
std::vector<float> vA_diphoton_reco_eta_;
std::vector<float> vA_diphoton_reco_phi_;

std::vector<float> vA_ditau_gen_m0_;
std::vector<float> vA_ditau_gen_dR_;
std::vector<float> vA_ditau_gen_E_;
std::vector<float> vA_ditau_gen_pT_;
std::vector<float> vA_ditau_gen_eta_;
std::vector<float> vA_ditau_gen_phi_;

// --- branches setup ---
void RecHitAnalyzer::branchesEvtSel_jet_h2aa2ditau_dipho(TTree* tree, edm::Service<TFileService>& fs) {
  tree->Branch("A_diphoton_gen_m0", &vA_diphoton_gen_m0_);
  tree->Branch("A_diphoton_gen_dR", &vA_diphoton_gen_dR_);
  tree->Branch("A_diphoton_gen_E", &vA_diphoton_gen_E_);
  tree->Branch("A_diphoton_gen_pT", &vA_diphoton_gen_pT_);
  tree->Branch("A_diphoton_gen_eta", &vA_diphoton_gen_eta_);
  tree->Branch("A_diphoton_gen_phi", &vA_diphoton_gen_phi_);
  tree->Branch("A_diphoton_reco_M", &vA_diphoton_reco_M_);
  tree->Branch("A_diphoton_reco_dR", &vA_diphoton_reco_dR_);
  tree->Branch("A_diphoton_reco_E", &vA_diphoton_reco_E_);
  tree->Branch("A_diphoton_reco_pT", &vA_diphoton_reco_pT_);
  tree->Branch("A_diphoton_reco_eta", &vA_diphoton_reco_eta_);
  tree->Branch("A_diphoton_reco_phi", &vA_diphoton_reco_phi_);

  tree->Branch("A_ditau_gen_m0", &vA_ditau_gen_m0_);
  tree->Branch("A_ditau_gen_dR", &vA_ditau_gen_dR_);
  tree->Branch("A_ditau_gen_E", &vA_ditau_gen_E_);
  tree->Branch("A_ditau_gen_pT", &vA_ditau_gen_pT_);
  tree->Branch("A_ditau_gen_eta", &vA_ditau_gen_eta_);
  tree->Branch("A_ditau_gen_phi", &vA_ditau_gen_phi_);
}

// --- event selection ---
bool RecHitAnalyzer::runEvtSel_jet_h2aa2ditau_dipho(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using GenPartColl = reco::GenParticleCollection;
  edm::Handle<GenPartColl> genParticles;
  iEvent.getByToken(genParticleCollectionT_, genParticles);

  edm::Handle<pat::PhotonCollection> photons;
  iEvent.getByToken(photonCollectionT_, photons);

  vAs_diphoton.clear();
  vAs_ditau.clear();
  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> vH(0,0,0,0);

  // find gen A->diphoton or A->ditau
  for (unsigned int iG = 0; iG < genParticles->size(); ++iG) {
      const auto& gen = genParticles->at(iG);
      int id = std::abs(gen.pdgId());
      if (id != 25 && id != 26) continue;
      if (gen.numberOfDaughters() != 2) continue;
      int dau0 = std::abs(gen.daughter(0)->pdgId());
      int dau1 = std::abs(gen.daughter(1)->pdgId());
      if (dau0 == 22 && dau1 == 22) {
          if (debug) std::cout << "Found A->#gamma#gamma idx=" << iG << " m=" << gen.mass() << std::endl;
          vAs_diphoton.push_back({iG, gen.pt()});
          vH += gen.p4();
      } else if (dau0 == 15 && dau1 == 15) {
          if (debug) std::cout << "Found A->#tau#tau idx=" << iG << " m=" << gen.mass() << std::endl;
          vAs_ditau.push_back({iG, gen.pt()});
          vH += gen.p4();
      }
  }

  // require exactly one of each
  if (vAs_diphoton.size() != 1 || vAs_ditau.size() != 1) return false;

  // record gen indices
  vGen_As_diphoton_Idxs.clear();
  vGen_As_ditau_Idxs.clear();
  vGen_As_diphoton_Idxs.push_back(vAs_diphoton[0].idx);
  vGen_As_ditau_Idxs.push_back(vAs_ditau[0].idx);

  // select reco photons near the ditau gen
  const auto& tauGen = genParticles->at(vAs_ditau[0].idx);
  double refEta = tauGen.eta(), refPhi = tauGen.phi();
  vPhotons.clear();
  for (unsigned int ip = 0; ip < photons->size(); ++ip) {
      const auto& ph = photons->at(ip);
      if (ph.pt() > pho_min_pT && reco::deltaR(ph.eta(), ph.phi(), refEta, refPhi) < 0.4) {
          vPhotons.push_back(ip);
      }
  }
  if (vPhotons.size() != 2) return false;

  vReco_First_Photons_Idxs.clear();
  vReco_Second_Photons_Idxs.clear();
  vReco_First_Photons_Idxs.push_back(vPhotons[0]);
  vReco_Second_Photons_Idxs.push_back(vPhotons[1]);

  return true;
}

// --- fill branches ---
void RecHitAnalyzer::fillEvtSel_jet_h2aa2ditau_dipho(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using GenPartColl = reco::GenParticleCollection;
  edm::Handle<GenPartColl> genParticles;
  iEvent.getByToken(genParticleCollectionT_, genParticles);

  edm::Handle<pat::PhotonCollection> photons;
  iEvent.getByToken(photonCollectionT_, photons);

  // clear all branch vectors
  vA_diphoton_gen_m0_.clear(); vA_diphoton_gen_dR_.clear(); vA_diphoton_gen_E_.clear();
  vA_diphoton_gen_pT_.clear(); vA_diphoton_gen_eta_.clear(); vA_diphoton_gen_phi_.clear();
  vA_diphoton_reco_M_.clear(); vA_diphoton_reco_dR_.clear(); vA_diphoton_reco_E_.clear();
  vA_diphoton_reco_pT_.clear(); vA_diphoton_reco_eta_.clear(); vA_diphoton_reco_phi_.clear();
  vA_ditau_gen_m0_.clear(); vA_ditau_gen_dR_.clear(); vA_ditau_gen_E_.clear();
  vA_ditau_gen_pT_.clear(); vA_ditau_gen_eta_.clear(); vA_ditau_gen_phi_.clear();

  // fill ditau gen
  for (auto idx : vGen_As_ditau_Idxs) {
      const auto& gen = genParticles->at(idx);
      vA_ditau_gen_E_.push_back(gen.energy());
      vA_ditau_gen_pT_.push_back(gen.pt());
      vA_ditau_gen_eta_.push_back(gen.eta());
      vA_ditau_gen_phi_.push_back(gen.phi());
      vA_ditau_gen_m0_.push_back(gen.mass());
      vA_ditau_gen_dR_.push_back(
          reco::deltaR(gen.daughter(0)->eta(), gen.daughter(0)->phi(),
                       gen.daughter(1)->eta(), gen.daughter(1)->phi()));
  }

  // fill diphoton gen
  for (auto idx : vGen_As_diphoton_Idxs) {
      const auto& gen = genParticles->at(idx);
      vA_diphoton_gen_E_.push_back(gen.energy());
      vA_diphoton_gen_pT_.push_back(gen.pt());
      vA_diphoton_gen_eta_.push_back(gen.eta());
      vA_diphoton_gen_phi_.push_back(gen.phi());
      vA_diphoton_gen_m0_.push_back(gen.mass());
      vA_diphoton_gen_dR_.push_back(
          reco::deltaR(gen.daughter(0)->eta(), gen.daughter(0)->phi(),
                       gen.daughter(1)->eta(), gen.daughter(1)->phi()));
  }

  // fill reconstructed diphoton
  for (size_t i = 0; i < vReco_First_Photons_Idxs.size(); ++i) {
      unsigned ip1 = vReco_First_Photons_Idxs[i];
      unsigned ip2 = vReco_Second_Photons_Idxs[i];
      const auto& ph1 = photons->at(ip1);
      const auto& ph2 = photons->at(ip2);
      auto diphoton = ph1.p4() + ph2.p4();
      vA_diphoton_reco_E_.push_back(diphoton.E());
      vA_diphoton_reco_pT_.push_back(diphoton.Pt());
      vA_diphoton_reco_eta_.push_back(diphoton.eta());
      vA_diphoton_reco_phi_.push_back(diphoton.phi());
      vA_diphoton_reco_M_.push_back(diphoton.M());
      vA_diphoton_reco_dR_.push_back(
          reco::deltaR(ph1.eta(), ph1.phi(), ph2.eta(), ph2.phi()));
  }
}