#include "MLAnalyzer_run3/RecHitAnalyzer/interface/RecHitAnalyzer.h"
#include <algorithm>
#include <tuple>

// -----------------------------------------------------------------------------
// Configuration
// -----------------------------------------------------------------------------
static constexpr unsigned int kMaxTrackTriplets = 8192; // keep <=5k
static constexpr float        kTripletPad       = -999.f;

// -----------------------------------------------------------------------------
// One std::vector<float> per channel & projection
//  * ECAL_tracksTriplet[proj]   : value (=1), eta, phi   -> ranked by track pT
//  * ECAL_tracksPtTriplet[proj] : pT,    eta, phi        -> ranked by pT
// size of each = 3*kMaxTrackTriplets
// -----------------------------------------------------------------------------
extern std::vector<float> vECAL_tracksTriplet_  [Nproj];
extern std::vector<float> vECAL_tracksPtTriplet_[Nproj];

#define DEFINE_TRIPLET_ARRAY(name) std::vector<float> name[Nproj];
DEFINE_TRIPLET_ARRAY(vECAL_tracksTriplet_)
DEFINE_TRIPLET_ARRAY(vECAL_tracksPtTriplet_)
#undef DEFINE_TRIPLET_ARRAY

namespace {
  struct TrackHit {
    float rank;   // what we sort by (pT)
    float val;    // stored as first element of triplet (1 or pT)
    float eta;
    float phi;
    bool operator<(const TrackHit& o) const { return rank > o.rank; } // descending
  };
}

void
RecHitAnalyzer::branchesTrackTripletsAtECAL (TTree *tree,
                                             edm::Service<TFileService>&)
{
  for (unsigned p=0; p<Nproj; ++p) {
    tree->Branch( (std::string("ECAL_tracks_triplet")   + projections[p]).c_str(),
                  &vECAL_tracksTriplet_[p]   );
    tree->Branch( (std::string("ECAL_tracksPt_triplet") + projections[p]).c_str(),
                  &vECAL_tracksPtTriplet_[p] );
  }
}

void
RecHitAnalyzer::fillTrackTripletsAtECAL (const edm::Event&  iEvent,
                                         const edm::EventSetup& iSetup,
                                         unsigned           proj)
{
  // Scratch buffers to collect every individual track after propagation
  std::vector<TrackHit> trackOcc;  // for ECAL_tracks     (val=1)
  std::vector<TrackHit> trackPt;   // for ECAL_tracksPt   (val=pT)

  // We repeat the propagation logic of the stitched module but
  // *do not* group into pixels – each hit kept separate.
  auto const& magfield    = iSetup.getData(magfieldToken_);
  auto const& caloGeom    = iSetup.getData(caloGeomToken_);
  auto const& transTrackB = iSetup.getData(transTrackBToken_);

  edm::Handle<reco::TrackCollection> tracksH;
  iEvent.getByToken(trackCollectionT_, tracksH);

  bool isPVgood=false;
  edm::Handle<reco::VertexCollection> pvColl;
  iEvent.getByToken(vertexCollectionT_, pvColl);
  isPVgood = pvColl.product()->size()>0;
  reco::Vertex thePV;
  if (isPVgood) thePV = pvColl.product()->at(0);

  reco::Track::TrackQuality tkQt = reco::Track::qualityByName("highPurity");

  for (const auto& tk : *tracksH) {

    if (tk.pt() <= 0.5)                continue;
    if (!tk.quality(tkQt))             continue;
    if (tk.charge() == 0)              continue;

    // Same "proj==4 => atPV only" selection as the stitched code
    if (proj==4) {
      bool pvMatch=false;
      for (const auto& tkPV : thePV.tracks()) {
        if ( std::abs(tk.pt()  - tkPV->pt())  < 1e-3 &&
             std::abs(tk.eta() - tkPV->eta()) < 1e-3 &&
             std::abs(tk.phi() - tkPV->phi()) < 1e-3 ) { pvMatch=true; break; }
      }
      if (!pvMatch) continue;
    }

    float eta  = 0.f;
    float phi  = 0.f;
    bool  ok   = false;

    switch (proj) {
      case 1: { eta = tk.eta(); phi = tk.phi(); ok=true; } break;
      case 2: case 0: case 4: case 5: {
        auto prop = spr::propagateTrackToECAL(&tk,&magfield);
        ok = prop.ok; if (ok){ eta=prop.direction.eta(); phi=prop.direction.phi(); }
      } break;
      case 3: {
        auto prop = spr::propagateTrackToHCAL(&tk,&magfield);
        ok = prop.ok; if (ok){ eta=prop.direction.eta(); phi=prop.direction.phi(); }
      } break;
    }

    if (!ok || std::abs(eta) > 3.) continue;
    // Require the trajectory to point to *any* ECAL crystal (barrel or endcap)
    DetId id( spr::findDetIdECAL(&caloGeom, eta, phi, false) );
    if (id == DetId(0)) continue; // safeguard

    // Rank order by pT
    TrackHit hOcc { static_cast<float>(tk.pt()), 1.f,
                static_cast<float>(eta),    static_cast<float>(phi) };

    TrackHit hPt  { static_cast<float>(tk.pt()), static_cast<float>(tk.pt()),
                static_cast<float>(eta),    static_cast<float>(phi) };

    trackOcc.emplace_back( std::move(hOcc) );
    trackPt.emplace_back ( std::move(hPt ) );
  }

  // ------------------------------------------------------------------
  // Sort & trim to <=5 000, then export to the fixed-length vectors
  // ------------------------------------------------------------------
  auto exportTriplets = [](std::vector<TrackHit>& src,
                           std::vector<float>&    dest)
  {
    std::sort(src.begin(), src.end());          // descending by rank
    const unsigned keep = std::min<unsigned>(src.size(), kMaxTrackTriplets);
    dest.assign(3*kMaxTrackTriplets, kTripletPad);

    for (unsigned i=0; i<keep; ++i) {
      dest[3*i]   = src[i].val;
      dest[3*i+1] = src[i].eta;
      dest[3*i+2] = src[i].phi;
    }
  };

  exportTriplets(trackOcc, vECAL_tracksTriplet_  [proj]);
  exportTriplets(trackPt , vECAL_tracksPtTriplet_[proj]);
}
