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

struct gen_obj {
  unsigned int idx;
  double pt;
};

// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesEvtSel_jet_h2aa2ditau_dipho ( TTree* tree, edm::Service<TFileService> &fs )
{
  tree->Branch("A_diphoton_gen_m0",         &vA_diphoton_gen_m0_);
  tree->Branch("A_diphoton_gen_dR",         &vA_diphoton_gen_dR_);
  tree->Branch("A_diphoton_gen_E",           &vA_diphoton_gen_E_);
  tree->Branch("A_diphoton_gen_pT",         &vA_diphoton_gen_pT_);
  tree->Branch("A_diphoton_gen_eta",       &vA_diphoton_gen_eta_);
  tree->Branch("A_diphoton_gen_phi",       &vA_diphoton_gen_phi_);
  tree->Branch("A_diphoton_reco_M",         &vA_diphoton_reco_M_);
  tree->Branch("A_diphoton_reco_dR",       &vA_diphoton_reco_dR_);
  tree->Branch("A_diphoton_reco_E",         &vA_diphoton_reco_E_);
  tree->Branch("A_diphoton_reco_pT",       &vA_diphoton_reco_pT_);
  tree->Branch("A_diphoton_reco_eta",     &vA_diphoton_reco_eta_);
  tree->Branch("A_diphoton_reco_phi",     &vA_diphoton_reco_phi_);

  tree->Branch("A_ditau_gen_m0",        &vA_ditau_gen_m0_);
  tree->Branch("A_ditau_gen_dR",        &vA_ditau_gen_dR_);
  tree->Branch("A_ditau_gen_E",          &vA_ditau_gen_E_);
  tree->Branch("A_ditau_gen_pT",        &vA_ditau_gen_pT_);
  tree->Branch("A_ditau_gen_eta",      &vA_ditau_gen_eta_);
  tree->Branch("A_ditau_gen_phi",      &vA_ditau_gen_phi_);
}

// Run event selection ___________________________________________________________________//
bool RecHitAnalyzer::runEvtSel_jet_h2aa2ditau_dipho ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleCollectionT_, genParticles);
  edm::Handle<reco::PhotonCollection> photons;
  iEvent.getByToken( photonCollectionT_, photons );

  vAs_diphoton.clear();
  vAs_ditau.clear();
  vPhotons.clear();
 
  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > vH;

  for ( unsigned int iG = 0; iG < genParticles->size(); iG++ ) {

    reco::GenParticleRef iGen( genParticles, iG );
    if ( std::abs(iGen->pdgId()) != 25 ||  std::abs(iGen->pdgId()) != 26) continue;
    if ( iGen->numberOfDaughters() != 2 ) continue;
    if ( abs(iGen->daughter(0)->pdgId()) == 22 || abs(iGen->daughter(1)->pdgId()) == 22 ) {
      if ( debug ) std::cout<<"*****************************************************"<< std::endl;
      if ( debug ) std::cout<< "iG:" << iG << " ID:" << iGen->pdgId() << " A->diphoton mass:" << iGen->mass() << std::endl;
      vAs_diphoton.push_back( iG );
      vH += iGen->p4();

    } else if ( abs(iGen->daughter(0)->pdgId()) == 15 || abs(iGen->daughter(1)->pdgId()) == 15 ) {
      if ( debug ) std::cout<<"*****************************************************"<< std::endl;
      if ( debug ) std::cout<< "iG:" << iG << " ID:" << iGen->pdgId() << " A->ditau mass:" << iGen->mass() << std::endl;
      gen_obj Gen_obj = { iG, std::abs(iGen->pt()) };
      vAs_ditau.push_back( iG );
      vH += iGen->p4();
    } else continue;

  }

  if ( vAs_diphoton.size() != 1 || vAs_ditau.size() != 1) return false;

  vGen_As_ditau_Idxs.clear();
  vGen_As_diphoton_Idxs.clear();
  reco::GenParticleRef iGen( genParticles, vAs_ditau[0] );
  if ( debug ) std::cout << " >> pT:" << iGen->pt() << " eta:" << iGen->eta() << " phi: " << iGen->phi() << " E:" << iGen->energy() << std::endl;

  vPhotons.clear()
  for (unsigned int iP = 0; iP < photons->size(); iP++) {
    reco::PhotonRef iP( photons, iP)
    if ( iP->pt() > pho_min_pT) && ( reco::deltaR( iP->eta(), iP->phi(), iGen->eta(), iGen->phi() ) < .4 ) {
      vPhotons.push_back( iP )
    }
  };
  
  if (vPhotons.size() != 2) return false;
  
  vReco_First_Photons_Idxs.push_back( vPhotons[0] );
  vReco_Second_Photons_Idxs.push_back( vPhotons[1] );


  vGen_As_ditau_Idxs.push_back( vAs_ditau[0] );
  reco::GenParticleRef iGen( genParticles, vAs_diphoton[0] );
  if ( debug ) std::cout << " >> pT:" << iGen->pt() << " eta:" << iGen->eta() << " phi: " << iGen->phi() << " E:" << iGen->energy() << std::endl;
  vGen_As_diphoton_Idxs.push_back( vAs_diphoton[0] );

  return true;
}

// Fill branches ___________________________________________________________________//
void RecHitAnalyzer::fillEvtSel_jet_h2aa2ditau_dipho ( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleCollectionT_, genParticles);

  std::Vector<reco::Photon> photons;
  iEvent.getByToken(photonCollectionT_, photons);

  //edm::Handle<PhotonCollection> photons;
  //iEvent.getByToken(photonCollectionT_, photons);

  vA_diphoton_gen_m0_.clear();
  vA_diphoton_gen_dR_.clear();
  vA_diphoton_gen_E_.clear();
  vA_diphoton_gen_pT_.clear();
  vA_diphoton_gen_eta_.clear();
  vA_diphoton_gen_phi_.clear();
  vA_diphoton_reco_M_.clear();
  vA_diphoton_reco_dR_.clear();
  vA_diphoton_reco_E_.clear();
  vA_diphoton_reco_pT_.clear();
  vA_diphoton_reco_eta_.clear();
  vA_diphoton_reco_phi_.clear();

  vA_ditau_gen_E_.clear();
  vA_ditau_gen_pT_.clear();
  vA_ditau_gen_eta_.clear();
  vA_ditau_gen_phi_.clear();
  vA_ditau_gen_m0_.clear();
  vA_ditau_gen_dR_.clear();

  vAs_diphoton.clear()
  vAs_ditau.clear()

  for ( unsigned int iG : vGen_As_ditau_Idxs ) {

    reco::GenParticleRef iGen( genParticles, iG );

    vA_ditau_gen_E_.push_back( std::abs(iGen->energy()) );
    vA_ditau_gen_pT_.push_back( std::abs(iGen->pt()) );
    vA_ditau_gen_eta_.push_back( iGen->eta() );
    vA_ditau_gen_phi_.push_back( iGen->phi() );
    vA_ditau_gen_m0_.push_back( iGen->mass() );
    vA_ditau_gen_dR_.push_back( reco::deltaR(iGen->daughter(0)->eta(),iGen->daughter(0)->phi(), iGen->daughter(1)->eta(),iGen->daughter(1)->phi()) );

  }

  for ( unsigned int iG : vGen_As_diphoton_Idxs ) {

    reco::GenParticleRef iGen( genParticles, iG );

    vA_diphoton_gen_E_.push_back( std::abs(iGen->energy()) );
    vA_diphoton_gen_pT_.push_back( std::abs(iGen->pt()) );
    vA_diphoton_gen_eta_.push_back( iGen->eta() );
    vA_diphoton_gen_phi_.push_back( iGen->phi() );
    vA_diphoton_gen_m0_.push_back( iGen->mass() );
    vA_diphoton_gen_dR_.push_back( reco::deltaR(iGen->daughter(0)->eta(),iGen->daughter(0)->phi(), iGen->daughter(1)->eta(),iGen->daughter(1)->phi()) );
  }

  for ( unsigned int iP = 0; iP < vReco_First_Photons_Idxs.size(); ++iP ) {
    unsigned int idx1 = vReco_First_Photons_Idxs[iP];
    unsigned int idx2 = vReco_Second_Photons_Idxs[iP];
  
    reco::PhotonRef photon1( photons, idx1 );
    reco::PhotonRef photon2( photons, idx2 );
  
    auto diphoton_A = photon1->p4() + photon2->p4();
  
    vA_diphoton_reco_E_ .push_back( std::abs(diphoton_A.energy()) );
    vA_diphoton_reco_pT_.push_back( std::abs(diphoton_A.pt())     );
    vA_diphoton_reco_eta_.push_back( diphoton_A.eta()             );
    vA_diphoton_reco_phi_.push_back( diphoton_A.phi()             );
    vA_diphoton_reco_M_ .push_back( diphoton_A.mass()            );
    vA_diphoton_reco_dR_.push_back(
      reco::deltaR(photon1->eta(), photon1->phi(),
                    photon2->eta(), photon2->phi())
    );
  }

}
