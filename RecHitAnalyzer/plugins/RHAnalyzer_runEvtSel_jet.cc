
#include "MLAnalyzer_run3/RecHitAnalyzer/interface/RecHitAnalyzer.h"

using std::vector;


// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesEvtSel_jet ( TTree* tree, edm::Service<TFileService> &fs ) {
  // Initialize event branches
  if ( task_ == "dijet_ditau" ) {
    branchesEvtSel_jet_h2aa2ditau_dipho( tree, fs );
  } 
}

// Run event selection ___________________________________________________________________//
bool RecHitAnalyzer::runEvtSel_jet ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {
  // Run event selection and reco-gen truth match
  if (task_ == "dijet_ditau") { 
    runEvtSel_jet_h2aa2ditau_dipho( iEvent, iSetup );
  }

  // Fill event branches
  if ( task_ == "dijet_ditau" ) {
    fillEvtSel_jet_h2aa2ditau_dipho( iEvent, iSetup );
  }

  return true;

}
