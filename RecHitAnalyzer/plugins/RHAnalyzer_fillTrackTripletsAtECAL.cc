#include "MLAnalyzer_run3/RecHitAnalyzer/interface/RecHitAnalyzer.h"
#include <algorithm>
#include <tuple>

// Configuration
static constexpr unsigned int kMaxTrackTriplets = 8192; // keep <=5k
static constexpr float        kTripletPad       = -999.f;
// Optional: make output eta/phi the exact crystal center
static constexpr bool         kSnapToCrystalCenter = true;

// One std::vector<float> per channel & projection
extern std::vector<float> vECAL_tracksTriplet_[Nproj];
extern std::vector<float> vECAL_tracksPtTriplet_[Nproj];

#define DEFINE_TRIPLET_ARRAY(name) std::vector<float> name[Nproj];
DEFINE_TRIPLET_ARRAY(vECAL_tracksTriplet_)
DEFINE_TRIPLET_ARRAY(vECAL_tracksPtTriplet_)
#undef DEFINE_TRIPLET_ARRAY

namespace {
  struct TrackHit {
    float     rank;     // sort key (pT)
    float     val;      // stored first (1 for occupancy, pT for pt-branch)
    float     eta;
    float     phi;
    uint32_t  ordinal;  // insertion order to break pT ties deterministically

    // Deterministic ordering: by pT desc, then by original insertion order
    bool operator<(const TrackHit& o) const {
      if (rank != o.rank) return rank > o.rank;   // descending pT
      return ordinal < o.ordinal;                 // stable tie-break
    }
  };

  inline float wrapPhi(float phi) {
    while (phi >  M_PI) phi -= 2.f*M_PI;
    while (phi <=-M_PI) phi += 2.f*M_PI;
    return phi;
  }
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
  auto padOut = [&](unsigned p) {
    vECAL_tracksTriplet_[p]   .assign(3*kMaxTrackTriplets, kTripletPad);
    vECAL_tracksPtTriplet_[p] .assign(3*kMaxTrackTriplets, kTripletPad);
  };

  // Scratch buffers
  std::vector<TrackHit> trackOcc;  // val=1
  std::vector<TrackHit> trackPt;   // val=pT

  // Conditions data
  auto const& magfield    = iSetup.getData(magfieldToken_);
  auto const& caloGeom    = iSetup.getData(caloGeomToken_);

  // Tracks
  edm::Handle<reco::TrackCollection> tracksH;
  iEvent.getByToken(trackCollectionT_, tracksH);
  if (!tracksH.isValid()) {
    padOut(proj);
    return;
  }
  trackOcc.reserve(tracksH->size());
  trackPt .reserve(tracksH->size());

  // PV (only needed for proj==4)
  bool isPVgood=false;
  edm::Handle<reco::VertexCollection> pvColl;
  iEvent.getByToken(vertexCollectionT_, pvColl);
  isPVgood = pvColl.isValid() && !pvColl->empty();
  reco::Vertex thePV;
  if (isPVgood) thePV = pvColl->at(0);

  if (proj==4 && !isPVgood) {
    padOut(proj);
    return;
  }

  const reco::Track::TrackQuality tkQt = reco::Track::qualityByName("highPurity");

  uint32_t ordinal = 0;
  for (const auto& tk : *tracksH) {
    const float pt  = tk.pt();
    if (pt <= 0.5f)                 continue;
    if (!tk.quality(tkQt))          continue;
    if (tk.charge() == 0)           continue;

    // proj==4 => require match to PV track list
    if (proj==4) {
      bool pvMatch=false;
      for (const auto& tkPV : thePV.tracks()) {
        if ( std::abs(tk.pt()  - tkPV->pt())  < 1e-3 &&
             std::abs(tk.eta() - tkPV->eta()) < 1e-3 &&
             std::abs(tk.phi() - tkPV->phi()) < 1e-3 ) { pvMatch=true; break; }
      }
      if (!pvMatch) continue;
    }

    // Determine (eta,phi) according to projection choice
    float eta  = 0.f;
    float phi  = 0.f;
    bool  ok   = false;

    // Require the trajectory to intersect an ECAL crystal at some point
    // for downstream alignment. Optionally snap to that crystal center
    float eta_for_id = 0.f, phi_for_id = 0.f;

    switch (proj) {
      case 1: {
        eta = tk.eta(); phi = tk.phi(); ok = true;
        auto prop = spr::propagateTrackToECAL(&tk, &magfield);
        if (!prop.ok) ok = false;
        else { eta_for_id = prop.direction.eta(); phi_for_id = prop.direction.phi(); }
      } break;
      case 2: case 0: case 4: case 5: {
        auto prop = spr::propagateTrackToECAL(&tk, &magfield);
        ok = prop.ok;
        if (ok) {
          eta = prop.direction.eta(); phi = prop.direction.phi();
          eta_for_id = eta; phi_for_id = phi;
        }
      } break;
      case 3: {
        auto prop = spr::propagateTrackToHCAL(&tk, &magfield);
        ok = prop.ok;
        if (ok) { eta = prop.direction.eta(); phi = prop.direction.phi(); }
        // For association check, also require an ECAL cell in front of HCAL:
        auto propE = spr::propagateTrackToECAL(&tk, &magfield);
        if (!(ok && propE.ok)) ok = false;
        else { eta_for_id = propE.direction.eta(); phi_for_id = propE.direction.phi(); }
      } break;
    }

    if (!ok) continue;
    if (std::abs(eta_for_id) > 3.f) continue;

    // Require a real ECAL crystal & optionally snap (eta,phi) to its center
    DetId id( spr::findDetIdECAL(&caloGeom, eta_for_id, phi_for_id, false) );
    if (id == DetId(0)) continue;

    if (kSnapToCrystalCenter) {
      auto const cell = caloGeom.getGeometry(id);  // shared_ptr<const CaloCellGeometry>
      if (!cell) continue;                         // check pointer is valid
      const GlobalPoint& gp = cell->getPosition();
      eta = gp.eta();
      phi = wrapPhi(gp.phi());
    } else {
      phi = wrapPhi(phi);
    }

    // Build paired entries with shared ordinal (ensures identical order post-sort)
    const uint32_t ord = ordinal++;
    trackOcc.emplace_back( TrackHit{ pt, 1.f,  eta, phi, ord } );
    trackPt .emplace_back( TrackHit{ pt, pt,  eta, phi, ord } );
  }

  // Sort both containers deterministically (same comparison), then export
  std::sort(trackOcc.begin(), trackOcc.end());
  std::sort(trackPt .begin(), trackPt .end());

  auto exportTriplets = [](const std::vector<TrackHit>& src,
                           std::vector<float>&          dest)
  {
    dest.assign(3*kMaxTrackTriplets, kTripletPad);
    const unsigned keep = std::min<unsigned>(src.size(), kMaxTrackTriplets);
    for (unsigned i=0; i<keep; ++i) {
      dest[3*i]   = src[i].val;  // 1 or pT depending on which vector we pass
      dest[3*i+1] = src[i].eta;
      dest[3*i+2] = src[i].phi;
    }
  };

  exportTriplets(trackOcc, vECAL_tracksTriplet_  [proj]);
  exportTriplets(trackPt , vECAL_tracksPtTriplet_[proj]);
}