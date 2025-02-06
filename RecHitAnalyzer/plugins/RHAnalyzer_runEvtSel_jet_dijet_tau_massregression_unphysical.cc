#include "MLAnalyzer_run3/RecHitAnalyzer/interface/RecHitAnalyzer.h"

using std::vector;

const unsigned nJets = 50; //TODO: use cfg level nJets_
TH1D *h_tau_mr_jet_pT_un;
TH1D *h_tau_mr_jet_E_un;
TH1D *h_tau_mr_jet_eta_un;
TH1D *h_tau_mr_jet_m0_un;
TH1D *h_tau_mr_jet_ma_un;
TH1D *h_tau_mr_jet_pta_un;
TH2D *h_tau_mr_jet_a_m_pt_un;
TH1D *h_tau_mr_jet_nJet_un;
TH1D *h_tau_mr_jet_isDiTau_un;
TH1D *h_tau_mr_jet_dR_un;
TH1D *h_tau_mr_jet_TaudR_un;
TH1D *h_tau_mr_jet_Tau1dR_un;
TH1D *h_tau_mr_jet_Tau2dR_un;
TH1D *h_tau_mr_jet_Tau1pT_un;
TH1D *h_tau_mr_jet_Tau2pT_un;
TH1D *h_tau_mr_jet_NrecoTaus_un;
TH1D *h_tau_mr_jet_NGenTaus_un;
TH1D *h_tau_mr_jet_recoTau1dR_un;
TH1D *h_tau_mr_jet_recoTau2dR_un;
TH1D *h_tau_mr_jet_n1dR_un;
TH1D *h_tau_mr_jet_n2dR_un;
vector<float> v_mr_jetIsDiTau_un;
vector<float> v_mr_jetadR_un;
vector<float> v_mr_ma_un;
vector<float> v_mr_pta_un;
vector<float> v_mr_jetTaudR_un;
vector<float> v_mr_jetTau1dR_un;
vector<float> v_mr_jetTau2dR_un;
vector<float> v_mr_jetTau1pT_un;
vector<float> v_mr_jetTau2pT_un;
vector<float> v_mr_jetNGenTaus_un;
vector<float> v_mr_jetNrecoTaus_un;
vector<float> v_mr_jetrecoTau1dR_un;
vector<float> v_mr_jetrecoTau2dR_un;
vector<float> v_mr_jetn1dR_un;
vector<float> v_mr_jetn2dR_un;

vector<float> v_mr_tau_jet_m0__un;
vector<float> v_mr_tau_jet_ma__un;
vector<float> v_mr_tau_jet_pta__un;
vector<float> v_mr_tau_jet_pt__un;
vector<float> v_mr_tau_jetPdgIds__un;
vector<float> v_mr_tau_jetIsDiTau__un;
vector<float> v_mr_tau_jetadR__un;
vector<float> v_mr_tau_jetTaudR__un;
vector<float> v_mr_tau_jetTau1dR__un;
vector<float> v_mr_tau_jetTau2dR__un;
vector<float> v_mr_tau_jetTau1pT__un;
vector<float> v_mr_tau_jetTau2pT__un;
vector<float> v_mr_tau_jetNGenTaus__un;
vector<float> v_mr_tau_jetNrecoTaus__un;
vector<float> v_mr_tau_jetrecoTau1dR__un;
vector<float> v_mr_tau_jetrecoTau2dR__un;
vector<float> v_mr_tau_jetn1dR__un;
vector<float> v_mr_tau_jetn2dR__un;

vector<float> v_mr_tau_subJetE__un[nJets];
vector<float> v_mr_tau_subJetPx__un[nJets];
vector<float> v_mr_tau_subJetPy__un[nJets];
vector<float> v_mr_tau_subJetPz__un[nJets];




// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesEvtSel_jet_dijet_tau_massregression_unphysical ( TTree* tree, edm::Service<TFileService> &fs ) {

  h_tau_mr_jet_E_un          = fs->make<TH1D>("h_jet_E"          , "E;E;Jets"                                   , 100,  0., 800.);
  h_tau_mr_jet_pT_un         = fs->make<TH1D>("h_jet_pT"         , "p_{T};p_{T};Jets"                           , 100,  0., 800.);
  h_tau_mr_jet_eta_un        = fs->make<TH1D>("h_jet_eta"        , "#eta;#eta;Jets"                             , 100, -5.,   5.);
  h_tau_mr_jet_nJet_un       = fs->make<TH1D>("h_jet_nJet"       , "nJet;nJet;Events"                           ,  10,  0.,  10.);
  h_tau_mr_jet_m0_un         = fs->make<TH1D>("h_jet_m0"         , "m_{jet};m_{jet};Jets"                       , 100,  0., 100.);
  h_tau_mr_jet_a_m_pt_un     = fs->make<TH2D>("h_a_m_pT"         , "m^{a} vs p_{T}^{a};m^{a} vs p_{T}^{a};Jets" ,  26, 0.,  3.6, 100, 30., 300.);
  h_tau_mr_jet_ma_un         = fs->make<TH1D>("h_jet_ma"         , "m^{a};m^{a};Jets"                           ,  26, 0.,  3.6);
  h_tau_mr_jet_pta_un        = fs->make<TH1D>("h_jet_pta"        , "p_{T}^{a};p_{T}^{a};Jets"                   ,  100, 0., 300.);
  h_tau_mr_jet_isDiTau_un    = fs->make<TH1D>("h_jet_isDiTau"    , "nIsDiTau;nIsDiTau;Jets"                     ,  10,  0.,  10.);
  h_tau_mr_jet_dR_un         = fs->make<TH1D>("h_jet_dR"         , "dR_{a,j};dR_{a,j};Jets"                     ,  50,  0.,  0.5);
  h_tau_mr_jet_TaudR_un      = fs->make<TH1D>("h_jet_TaudR"      , "dR_{#tau,#tau};dR_{#tau,#tau};Jets"         ,  50,  0.,   1.);
  h_tau_mr_jet_Tau1dR_un     = fs->make<TH1D>("h_jet_Tau1dR"     , "dR_{#tau_{1},j};dR_{#tau_{1},j};Jets"       ,  50,  0.,  0.5);
  h_tau_mr_jet_Tau2dR_un     = fs->make<TH1D>("h_jet_Tau2dR"     , "dR_{#tau_{2},j};dR_{#tau_{2},j};Jets"       ,  50,  0.,  0.5);
  h_tau_mr_jet_Tau1pT_un     = fs->make<TH1D>("h_jet_Tau1pT"     , "p_{T}^{#tau_{1}};p_{T}^{#tau_{1}};Jets"     ,  50,  0.,  100);
  h_tau_mr_jet_Tau2pT_un     = fs->make<TH1D>("h_jet_Tau2pT"     , "p_{T}^{#tau_{2}};p_{T}^{#tau_{2}};Jets"     ,  50,  0.,  100);
  h_tau_mr_jet_NGenTaus_un  = fs->make<TH1D>("h_jet_NGenTaus"    , "N#tau^{RECO};N#tau^{RECO};Jets"             ,   5,  0.,   5.);
  h_tau_mr_jet_NrecoTaus_un  = fs->make<TH1D>("h_jet_NrecoTaus"  , "N#tau^{RECO};N#tau^{RECO};Jets"             ,   5,  0.,   5.);
  h_tau_mr_jet_recoTau1dR_un = fs->make<TH1D>("h_jet_recoTau1dR" , "dR_{#tau_{1}^{RECO},j};dR_{#tau_{1}^{RECO},j};Jets" ,  50,  0.,  0.5);
  h_tau_mr_jet_recoTau2dR_un = fs->make<TH1D>("h_jet_recoTau2dR" , "dR_{#tau_{2}^{RECO},j};dR_{#tau_{2}^{RECO},j};Jets" ,  25,  0.,  0.5);
  h_tau_mr_jet_n1dR_un       = fs->make<TH1D>("h_jet_n1dR"       , "dR_{#eta_{1},j};dR_{#eta_{1},j};Jets"       ,  25,  0.,  0.5);
  h_tau_mr_jet_n2dR_un       = fs->make<TH1D>("h_jet_n2dR"       , "dR_{#eta_{2},j};dR_{#eta_{2},j};Jets"       ,  25,  0.,  0.5);

  tree->Branch("jetM_un",       &v_mr_tau_jet_m0__un);
  tree->Branch("jetPt_un",      &v_mr_tau_jet_pt__un);
  tree->Branch("jetPdgIds_un",  &v_mr_tau_jetPdgIds__un);
  tree->Branch("jetadR_un",     &v_mr_tau_jetadR__un);
  tree->Branch("jetIsDiTau_un", &v_mr_tau_jetIsDiTau__un);
  tree->Branch("a_m_un",        &v_mr_tau_jet_ma__un);
  tree->Branch("a_pt_un",       &v_mr_tau_jet_pta__un);
  tree->Branch("jetpT_un",      &v_mr_tau_jet_pt__un);
  tree->Branch("TaudR_un",      &v_mr_tau_jetTaudR__un);
  tree->Branch("Tau1dR_un",     &v_mr_tau_jetTau1dR__un);
  tree->Branch("Tau2dR_un",     &v_mr_tau_jetTau2dR__un);
  tree->Branch("Tau1pT_un",     &v_mr_tau_jetTau1pT__un);
  tree->Branch("Tau2pT_un",     &v_mr_tau_jetTau2pT__un);
  tree->Branch("NGenTaus_un",   &v_mr_tau_jetNGenTaus__un);
  tree->Branch("NrecoTaus_un",  &v_mr_tau_jetNrecoTaus__un);
  tree->Branch("recoTau1dR_un", &v_mr_tau_jetrecoTau1dR__un);
  tree->Branch("recoTau2dR_un", &v_mr_tau_jetrecoTau2dR__un);
  tree->Branch("n1dR_un",       &v_mr_tau_jetn1dR__un);
  tree->Branch("n2dR_un",       &v_mr_tau_jetn2dR__un);

  char hname[50];
  for ( unsigned iJ = 0; iJ != nJets; iJ++ ) {
    sprintf(hname, "subJet%d_E", iJ);
    tree->Branch(hname,            &v_mr_tau_subJetE__un[iJ]);
    sprintf(hname, "subJet%d_Px", iJ);
    tree->Branch(hname,            &v_mr_tau_subJetPx__un[iJ]);
    sprintf(hname, "subJet%d_Py", iJ);
    tree->Branch(hname,            &v_mr_tau_subJetPy__un[iJ]);
    sprintf(hname, "subJet%d_Pz", iJ);
    tree->Branch(hname,            &v_mr_tau_subJetPz__un[iJ]);
  }

} // branchesEvtSel_jet_dijet_tau_massregression()

// Run jet selection _____________________________________________________//
bool RecHitAnalyzer::runEvtSel_jet_dijet_tau_massregression_unphysical( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken( genParticleCollectionT_, genParticles );

  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);

  edm::Handle<reco::PFTauCollection> taus;
  iEvent.getByToken(tauCollectionT_, taus);

  vJetIdxs.clear();
  v_mr_tau_jetPdgIds__un.clear();
  v_mr_jetIsDiTau_un.clear();
  v_mr_jetadR_un.clear();
  v_mr_ma_un.clear();
  v_mr_pta_un.clear();
  v_mr_jetTaudR_un.clear();
  v_mr_jetTau1dR_un.clear();
  v_mr_jetTau2dR_un.clear();
  v_mr_jetTau1pT_un.clear();
  v_mr_jetTau2pT_un.clear();
  v_mr_jetNGenTaus_un.clear();
  v_mr_jetNrecoTaus_un.clear();
  v_mr_jetrecoTau1dR_un.clear();
  v_mr_jetrecoTau2dR_un.clear();
  v_mr_jetn1dR_un.clear();
  v_mr_jetn2dR_un.clear();

  /*
  edm::Handle<reco::GenJetCollection> genJets;
  iEvent.getByToken(genJetCollectionT_, genJets);
  std::vector<float> v_mr_tau_jetFakePhoIdxs;
  */

  unsigned int nMatchedJets_un = 0;
  unsigned int nMatchedRecoTaus_un = 0;
  unsigned int aPdgId_un           = 0;
  bool MatchedPseudoScalar_un = false;
  float a_mass_un = -99.;
  float a_pt_un   = -99.;
  float dRa_un    = -99.;
  float tausdR_un =  99.;
  float tau1dR_un = -99.;
  float tau2dR_un = -99.;
  float tau1pT_un = -99.;
  float tau2pT_un = -99.;
  float recotau1dR_un = -99.;
  float recotau2dR_un = -99.;
  float n1dR_un   = -99.;
  float n2dR_un   = -99.;



  if ( debug ) std::cout << " >>>>>>>>>>>>>>>>>>>> evt:" << std::endl;
  if ( debug ) std::cout << " JETS IN THE EVENT = " << jets->size() << " | Selection requires minpT = " << minJetPt_ << " and maxEta = "<< maxJetEta_ << std::endl;
  // Loop over jets
  for ( unsigned iJ(0); iJ != jets->size(); ++iJ ) {
    reco::PFJetRef iJet( jets, iJ );
    if (debug ) std::cout << "  >>>>>> Jet [" << iJ << "] -> Pt: " << iJet->pt() << ", Eta: " << iJet->eta() << ", Phi: " << iJet->phi() << std::endl;
    if ( std::abs(iJet->pt())  < minJetPt_ ) continue;
    if ( std::abs(iJet->eta()) > maxJetEta_ ) continue;
    unsigned int nMatchedGenParticles_un = 0;
    bool passedGenSel_un = false;
    unsigned int iGenParticle_un = 0;
    unsigned int NTau1Daughters_un = 0;
    unsigned int NTau2Daughters_un = 0;
    for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen) {
      // if ( abs(iGen->pdgId()) != 25 || abs(iGen->daughter(0)->pdgId()) != 15 || abs(iGen->daughter(1)->pdgId()) != 16 ) continue;
      if ( abs(iGen->pdgId()) != 15 || iGen->numberOfDaughters() < 2 ) continue;

      if (iGen->daughter(0)->numberOfMothers() != 1 || iGen->daughter(1)->numberOfMothers() != 1) continue;
      ++iGenParticle_un;
      // float dR_un = reco::deltaR( iJet->eta(),iJet->phi(), iGen->daughter(0)->eta(),iGen->daughter(0)->phi() );
      float dR_un = reco::deltaR( iJet->eta(),iJet->phi(), iGen->eta(),iGen->phi() );
      if ( dR_un > 0.4 ) continue;

      // if ( iGen->numberOfMothers() != 1 ) continue;
      aPdgId_un = std::abs(iGen->pdgId());
      // if ( abs(iGen->pdgId()) == 25 && iGen->mass() < 15 && iGen->mass() > 3.5) {
      if (iGen->mass() < 3.601 && iGen->mass() > 1.799) {
        MatchedPseudoScalar_un = true;
      }
      else continue;

      if (nMatchedGenParticles_un == 0) {

        a_mass_un = iGen->mass();
        a_pt_un   = iGen->pt();
        dRa_un = reco::deltaR( iJet->eta(),iJet->phi(), iGen->eta(),iGen->phi() );
        // if ( iGen->numberOfDaughters() == 2 ){
        if ( iGen->numberOfDaughters() > 1 ){
          // if (abs(iGen->daughter(0)->pdgId()) == 15 && abs(iGen->daughter(1)->pdgId()) == 16){
          if (abs(iGen->daughter(0)->pdgId()) > 15 && abs(iGen->daughter(1)->pdgId()) > 15){
            tausdR_un = reco::deltaR( iGen->daughter(0)->eta(),iGen->daughter(0)->phi(), iGen->daughter(1)->eta(),iGen->daughter(1)->phi() );
            // if ( tausdR > 0.4 ) continue;

            if ( debug ) std::cout << "   TAUS   -> status: " << iGen->status() << ", id: " << iGen->pdgId() << ", nDaught: " << iGen->numberOfDaughters() << " | pt: "<< iGen->pt() << " eta: " <<iGen->eta() << " phi: " <<iGen->phi() << " mass: " <<iGen->mass() << std::endl;



            if ( iGen->daughter(0)->pt() > iGen->daughter(1)->pt() ) {
              tau1pT_un = iGen->daughter(0)->pt();
              tau2pT_un = iGen->daughter(1)->pt();
              tau1dR_un = reco::deltaR( iJet->eta(),iJet->phi(), iGen->daughter(0)->eta(),iGen->daughter(0)->phi() );
              tau2dR_un = reco::deltaR( iJet->eta(),iJet->phi(), iGen->daughter(1)->eta(),iGen->daughter(1)->phi() );

              if ( debug ) std::cout << "   >>>>>> Taus dR = " << tausdR_un << " , tau1 dR = " << tau1dR_un << " , tau2 dR = " << tau2dR_un << std::endl;
              // NTau1Daughters_un = iGen->daughter(0)->numberOfDaughters();
              // NTau2Daughters_un = iGen->daughter(1)->numberOfDaughters();
              NTau1Daughters_un = iGen->numberOfDaughters();
              NTau2Daughters_un = iGen->numberOfDaughters();
              if ( debug ) std::cout << "    >>>>>> # Tau 1 daughters = " << NTau1Daughters_un << ",  # Tau 2 daughters = "<< NTau2Daughters_un << std::endl;
              for (unsigned int iDaughter_un = 0; iDaughter_un != NTau1Daughters_un; ++iDaughter_un){
                if ( debug ) std::cout << "   >>>>>> Tau 1 Decay = " << iGen->daughter(iDaughter_un)->pdgId() << std::endl;
                if ( abs(iGen->daughter(iDaughter_un)->pdgId()) == 16 ) n1dR_un = reco::deltaR( iJet->eta(),iJet->phi(), iGen->daughter(iDaughter_un)->eta(), iGen->daughter(iDaughter_un)->phi() );
              }
              // for (unsigned int iDaughter_un = 0; iDaughter_un != NTau2Daughters_un; ++iDaughter_un){
              //   if ( debug ) std::cout << "   >>>>>> Tau 2 Decay = " << iGen->daughter(iDaughter_un)->pdgId() << std::endl;
              //   if ( abs(iGen->daughter(iDaughter_un)->pdgId()) == 16 ) n2dR_un = reco::deltaR( iJet->eta(),iJet->phi(), iGen->daughter(1)->daughter(iDaughter_un)->eta(), iGen->daughter(1)->daughter(iDaughter_un)->phi() );
              // }
            } else {
              tau1pT_un = iGen->daughter(1)->pt();
              tau2pT_un = iGen->daughter(0)->pt();
              tau1dR_un = reco::deltaR( iJet->eta(),iJet->phi(), iGen->daughter(1)->eta(),iGen->daughter(1)->phi() );
              tau2dR_un = reco::deltaR( iJet->eta(),iJet->phi(), iGen->daughter(0)->eta(),iGen->daughter(0)->phi() );

              if ( debug ) std::cout << "   >>>>>> Taus dR = " << tausdR_un << " , tau1 dR = " << tau1dR_un << " , tau2 dR = " << tau2dR_un << std::endl;
              NTau1Daughters_un = iGen->numberOfDaughters();
              NTau2Daughters_un = iGen->numberOfDaughters();
              if ( debug ) std::cout << "    >>>>>> # Tau 1 daughters = " << NTau1Daughters_un << ",  # Tau 2 daughters = "<< NTau2Daughters_un << std::endl;
              for (unsigned int iDaughter_un = 0; iDaughter_un != NTau1Daughters_un; ++iDaughter_un){
                if ( debug ) std::cout << "   >>>>>> Tau 1 Decay = " << iGen->daughter(iDaughter_un)->pdgId() << std::endl;
                if ( abs(iGen->daughter(iDaughter_un)->pdgId()) == 16 ) n1dR_un = reco::deltaR( iJet->eta(),iJet->phi(), iGen->daughter(iDaughter_un)->eta(), iGen->daughter(iDaughter_un)->phi() );
              }
              // for (unsigned int iDaughter_un = 0; iDaughter_un != NTau2Daughters_un; ++iDaughter_un){
              //   if ( debug ) std::cout << "   >>>>>> Tau 2 Decay = " << iGen->daughter(0)->daughter(iDaughter_un)->pdgId() << std::endl;
              //   if ( abs(iGen->daughter(0)->daughter(iDaughter_un)->pdgId()) == 16 ) n2dR_un = reco::deltaR( iJet->eta(),iJet->phi(), iGen->daughter(0)->daughter(iDaughter_un)->eta(), iGen->daughter(0)->daughter(iDaughter_un)->phi() );
              // }
            } // end else pt2 > pt1
          }
        }
        //}
        if (debug ) std::cout << "   >>>>>> n1 dR = " << n1dR_un << " , n2 dR = " << n2dR_un << std::endl;
      }

      if ( debug ) std::cout << "   GEN particle " << iGenParticle_un << " -> status: " << iGen->status() << ", id: " << iGen->pdgId() << ", nDaught: " << iGen->numberOfDaughters() << " nMoms: " <<iGen->numberOfMothers() << " | pt: "<< iGen->pt() << " eta: " <<iGen->eta() << " phi: " <<iGen->phi() << " | dR = "<< dR_un << std::endl ;

      if (debug ) std::cout << "  >>>>>> Jet [" << iJ << "] matched particle [" << iGenParticle_un << "] -> pdgId: " << std::abs(iGen->pdgId()) << " | dR: " << dR_un << "| Pt: " << iJet->pt() << ", Eta: " << iJet->eta() << ", Phi: " << iJet->phi() << std::endl;

      ++nMatchedGenParticles_un;
    } // primary gen particles
    if ( nMatchedGenParticles_un > 0 ) passedGenSel_un = true;
    if (passedGenSel_un) {
      ++nMatchedJets_un;

      //Lookin at RecoTaus
      if (taus->size() == 0) {
        if (debug ) std::cout << "   !!!!!!!!!!  NO RECO TAUS IN THIS EVENT  !!!!!!!!!!"<< std::endl;
      }

      for ( unsigned iT(0); iT != taus->size(); ++iT ) {
        reco::PFTauRef iTau( taus, iT );
        float recotaudR_un = reco::deltaR( iJet->eta(),iJet->phi(), iTau->eta(),iTau->phi() );
        if ( recotaudR_un < 0.4 && nMatchedRecoTaus_un == 0 ) {
          if ( debug ) std::cout << "Reco Tau [" << iT << "]  matched jet [" << iJ << "] : dR = " << recotaudR_un << std::endl;
          recotau1dR_un = recotaudR_un;
          ++nMatchedRecoTaus_un;
        } else if ( recotaudR_un < 0.4 && nMatchedRecoTaus_un == 1 ) {
          if ( debug ) std::cout << "Reco Tau [" << iT << "]  matched jet [" << iJ << "] : dR = " << recotaudR_un << std::endl;
          if (recotaudR_un < recotau1dR_un) {
            recotau2dR_un = recotau1dR_un;
            recotau1dR_un = recotaudR_un;
          } else recotau2dR_un = recotaudR_un;
          ++nMatchedRecoTaus_un;
        } else if ( debug && recotaudR_un < 0.4 && nMatchedRecoTaus_un > 1 ) {
          std::cout << "   !!!!!!!!!!  FOUND MORE THAN 2 TAUS INSIDE JET CONE OF 0.4 !!!!!!!!!!"<< std::endl;
          if (recotaudR_un < recotau2dR_un && recotaudR_un < recotau1dR_un) {
            if (recotau1dR_un < recotau2dR_un) recotau2dR_un = recotau1dR_un;
            recotau1dR_un = recotaudR_un;
          } else if (recotaudR_un < recotau2dR_un && recotaudR_un > recotau1dR_un) recotau2dR_un = recotaudR_un;
          ++nMatchedRecoTaus_un;
        } else if ( debug ) {
          std::cout << "   !!!!!!!!!!  NO MATCH FOR Reco Tau [" << iT << "]  with jet [" << iJ << "] : dR = " << recotaudR_un << std::endl;
        }
      }

      vJetIdxs.push_back( iJ );
      v_mr_tau_jetPdgIds__un.push_back( aPdgId_un );
      v_mr_jetadR_un.push_back( dRa_un );
      v_mr_jetTaudR_un.push_back( tausdR_un );
      v_mr_ma_un.push_back( a_mass_un );
      v_mr_pta_un.push_back( a_pt_un );
      v_mr_jetTau1dR_un.push_back( tau1dR_un );
      v_mr_jetTau2dR_un.push_back( tau2dR_un );
      v_mr_jetTau1pT_un.push_back( tau1pT_un );
      v_mr_jetTau2pT_un.push_back( tau2pT_un );
      v_mr_jetn1dR_un.push_back( n1dR_un );
      v_mr_jetn2dR_un.push_back( n2dR_un );
      v_mr_jetNGenTaus_un.push_back( nMatchedGenParticles_un );
      v_mr_jetNrecoTaus_un.push_back( nMatchedRecoTaus_un );
      v_mr_jetrecoTau1dR_un.push_back( recotau1dR_un );
      v_mr_jetrecoTau2dR_un.push_back( recotau2dR_un );
      v_mr_jetIsDiTau_un.push_back( MatchedPseudoScalar_un );

    }

  } // reco jets
  if ( debug ) std::cout << " Matched reco jets " << nMatchedJets_un << std::endl;
  if ( debug ) std::cout << " Matched reco Taus " << nMatchedRecoTaus_un << std::endl;

  // Check jet multiplicity
  if ( nMatchedJets_un < 1 && nMatchedRecoTaus_un < 1) return false;  // modified such that events has at lease one reco jet match to gen tau and reco tau

  if ( debug ) std::cout << " >> has_jet_dijet_tau_massregression_unphysical: passed" << std::endl;
  return true;

} // runEvtSel_jet_dijet_tau_massregression()

// Fill branches and histograms _____________________________________________________//
void RecHitAnalyzer::fillEvtSel_jet_dijet_tau_massregression_unphysical ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);

  h_tau_mr_jet_nJet_un->Fill( vJetIdxs.size() );

  v_mr_tau_jet_pt__un.clear();
  v_mr_tau_jet_m0__un.clear();
  v_mr_tau_jet_ma__un.clear();
  v_mr_tau_jet_pta__un.clear();
  v_mr_tau_jetIsDiTau__un.clear();
  v_mr_tau_jetadR__un.clear();
  v_mr_tau_jetTaudR__un.clear();
  v_mr_tau_jetTau1dR__un.clear();
  v_mr_tau_jetTau2dR__un.clear();
  v_mr_tau_jetTau1pT__un.clear();
  v_mr_tau_jetTau2pT__un.clear();
  v_mr_tau_jetNGenTaus__un.clear();
  v_mr_tau_jetNrecoTaus__un.clear();
  v_mr_tau_jetrecoTau1dR__un.clear();
  v_mr_tau_jetrecoTau2dR__un.clear();
  v_mr_tau_jetn1dR__un.clear();
  v_mr_tau_jetn2dR__un.clear();

  for ( unsigned iJ(0); iJ != vJetIdxs.size(); ++iJ ) {

    reco::PFJetRef iJet( jets, vJetIdxs[iJ] );

    // Fill histograms
    h_tau_mr_jet_pT_un->Fill( std::abs(iJet->pt()) );
    h_tau_mr_jet_eta_un->Fill( iJet->eta() );
    h_tau_mr_jet_E_un->Fill( iJet->energy() );
    h_tau_mr_jet_m0_un->Fill( iJet->mass() );
    h_tau_mr_jet_isDiTau_un->Fill( v_mr_jetIsDiTau_un[iJ] );
    h_tau_mr_jet_dR_un->Fill( v_mr_jetadR_un[iJ] );
    h_tau_mr_jet_ma_un->Fill( v_mr_ma_un[iJ] );
    h_tau_mr_jet_pta_un->Fill( v_mr_pta_un[iJ] );
    h_tau_mr_jet_a_m_pt_un->Fill( v_mr_ma_un[iJ], v_mr_pta_un[iJ] );
    h_tau_mr_jet_TaudR_un->Fill( v_mr_jetTaudR_un[iJ] );
    h_tau_mr_jet_Tau1dR_un->Fill( v_mr_jetTau1dR_un[iJ] );
    h_tau_mr_jet_Tau2dR_un->Fill( v_mr_jetTau2dR_un[iJ] );
    h_tau_mr_jet_Tau1pT_un->Fill( v_mr_jetTau1pT_un[iJ] );
    h_tau_mr_jet_Tau2pT_un->Fill( v_mr_jetTau2pT_un[iJ] );
    h_tau_mr_jet_NGenTaus_un->Fill( v_mr_jetNGenTaus_un[iJ] );
    h_tau_mr_jet_NrecoTaus_un->Fill( v_mr_jetNrecoTaus_un[iJ] );
    h_tau_mr_jet_recoTau1dR_un->Fill( v_mr_jetrecoTau1dR_un[iJ] );
    h_tau_mr_jet_recoTau2dR_un->Fill( v_mr_jetrecoTau2dR_un[iJ] );
    h_tau_mr_jet_n1dR_un->Fill( v_mr_jetn1dR_un[iJ] );
    h_tau_mr_jet_n2dR_un->Fill( v_mr_jetn2dR_un[iJ] );

    // Fill branches
    v_mr_tau_jet_pt__un.push_back( iJet->pt() );
    v_mr_tau_jet_m0__un.push_back( iJet->mass() );
    v_mr_tau_jet_ma__un.push_back( v_mr_ma_un[iJ] );
    v_mr_tau_jet_pta__un.push_back( v_mr_pta_un[iJ] );
    v_mr_tau_jetIsDiTau__un.push_back( v_mr_jetIsDiTau_un[iJ] );
    v_mr_tau_jetadR__un.push_back( v_mr_jetadR_un[iJ] );
    v_mr_tau_jetTaudR__un.push_back( v_mr_jetTaudR_un[iJ] );
    v_mr_tau_jetTau1dR__un.push_back( v_mr_jetTau1dR_un[iJ] );
    v_mr_tau_jetTau2dR__un.push_back( v_mr_jetTau2dR_un[iJ] );
    v_mr_tau_jetTau1pT__un.push_back( v_mr_jetTau1pT_un[iJ] );
    v_mr_tau_jetTau2pT__un.push_back( v_mr_jetTau2pT_un[iJ] );
    v_mr_tau_jetNGenTaus__un.push_back( v_mr_jetNGenTaus_un[iJ] );
    v_mr_tau_jetNrecoTaus__un.push_back( v_mr_jetNrecoTaus_un[iJ] );
    v_mr_tau_jetrecoTau1dR__un.push_back( v_mr_jetrecoTau1dR_un[iJ] );
    v_mr_tau_jetrecoTau2dR__un.push_back( v_mr_jetrecoTau2dR_un[iJ] );
    v_mr_tau_jetn1dR__un.push_back( v_mr_jetn1dR_un[iJ] );
    v_mr_tau_jetn2dR__un.push_back( v_mr_jetn2dR_un[iJ] );

    // Gen jet constituents
    v_mr_tau_subJetE__un[iJ].clear();
    v_mr_tau_subJetPx__un[iJ].clear();
    v_mr_tau_subJetPy__un[iJ].clear();
    v_mr_tau_subJetPz__un[iJ].clear();
    //std::vector<reco::PFCandidatePtr> jetConstituents = iJet->getPFConstituents();
    unsigned int nConstituents_un = iJet->getPFConstituents().size();
    if ( debug ) std::cout << " >>>>>> Jet " << iJ << " has "<< nConstituents_un << " constituents" << std::endl;
    for ( unsigned int j = 0; j < nConstituents_un; j++ ) {
      const reco::PFCandidatePtr subJet = iJet->getPFConstituent( j );
      if ( debug ) std::cout << " >>>>>>>>>>>  Jet constituent " << j << "-> E:" << subJet->energy() << " px:" << subJet->px() << " py:" << subJet->py() << " pz:" << subJet->pz() << std::endl;
      v_mr_tau_subJetE__un[iJ].push_back( subJet->energy() );
      v_mr_tau_subJetPx__un[iJ].push_back( subJet->px() );
      v_mr_tau_subJetPy__un[iJ].push_back( subJet->py() );
      v_mr_tau_subJetPz__un[iJ].push_back( subJet->pz() );
    }
  }

} // fillEvtSel_jet_dijet_tau_massregression_unphysica()
