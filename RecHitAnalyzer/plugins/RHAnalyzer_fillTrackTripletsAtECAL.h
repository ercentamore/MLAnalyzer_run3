#ifndef RHAnalyzer_fillTrackTripletsAtECAL_h
#define RHAnalyzer_fillTrackTripletsAtECAL_h 1

#include "MLAnalyzer_run3/RecHitAnalyzer/interface/RecHitAnalyzer.h"

// -----------------------------------------------------------------------------
// Configuration
// -----------------------------------------------------------------------------
static constexpr unsigned int kMaxTrackTriplets = 5000; // keep <=5k
static constexpr float        kTripletPad       = -999.f;

// -----------------------------------------------------------------------------
// One std::vector<float> per channel & projection
//  * ECAL_tracksTriplet[proj]   : value (=1), eta, phi   -> ranked by track pT
//  * ECAL_tracksPtTriplet[proj] : pT,    eta, phi        -> ranked by pT
// size of each = 3*kMaxTrackTriplets
// -----------------------------------------------------------------------------
extern std::vector<float> vECAL_tracksTriplet_  [Nproj];
extern std::vector<float> vECAL_tracksPtTriplet_[Nproj];

// -----------------------------------------------------------------------------
// Function declarations
// -----------------------------------------------------------------------------
void RecHitAnalyzer::branchesTrackTripletsAtECAL   (TTree*, edm::Service<TFileService>&);
void RecHitAnalyzer::fillTrackTripletsAtECAL       (const edm::Event&,
                                                    const edm::EventSetup&,
                                                    unsigned int proj);

#endif
