#ifndef RHAnalyzer_fillTRKTriplets_h
#define RHAnalyzer_fillTRKTriplets_h 1

#include "MLAnalyzer_run3/RecHitAnalyzer/interface/RecHitAnalyzer.h"

// -----------------------------------------------------------------------------
// Constants
// -----------------------------------------------------------------------------
static constexpr unsigned int kMaxTripletHits = 5000;   // keep at most 5k
static constexpr float        kPadValue       = -999.f; // pad sentinel

// -----------------------------------------------------------------------------
// Per‑detector containers  (value, eta, phi packed as V,E,P,V,E,P, ...)
// size will always be 3*kMaxTripletHits
// -----------------------------------------------------------------------------
extern std::vector<float> vBPIX_TRKTriplets_[nBPIX][Nhitproj];
extern std::vector<float> vFPIX_TRKTriplets_[nFPIX][Nhitproj];
extern std::vector<float> vTIB_TRKTriplets_ [nTIB ][Nhitproj];
extern std::vector<float> vTID_TRKTriplets_ [nTID ][Nhitproj];
extern std::vector<float> vTOB_TRKTriplets_ [nTOB ][Nhitproj];
extern std::vector<float> vTEC_TRKTriplets_ [nTEC ][Nhitproj];

// -----------------------------------------------------------------------------
// Function declarations
// -----------------------------------------------------------------------------
void RecHitAnalyzer::branchesTRKTriplets      (TTree            *tree,
                                                     edm::Service<TFileService>&);
void RecHitAnalyzer::fillTRKTriplets          (const edm::Event&,
                                                     const edm::EventSetup&,
                                                     unsigned int proj);

#endif
