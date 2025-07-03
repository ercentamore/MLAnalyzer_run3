#include "RHAnalyzer_fillTRKTriplets.h"
#include <algorithm>   // std::sort
#include <tuple>

// -----------------------------------------------------------------------------
// Container definitions (one vector / layer / projection)
// -----------------------------------------------------------------------------
#define DEFINE_TRIPLET_CONTAINER(name, NLAYER)         \
  std::vector<float> name[NLAYER][Nhitproj];

DEFINE_TRIPLET_CONTAINER(vBPIX_TRKTriplets_, nBPIX)
DEFINE_TRIPLET_CONTAINER(vFPIX_TRKTriplets_, nFPIX)
DEFINE_TRIPLET_CONTAINER(vTIB_TRKTriplets_,  nTIB )
DEFINE_TRIPLET_CONTAINER(vTID_TRKTriplets_,  nTID )
DEFINE_TRIPLET_CONTAINER(vTOB_TRKTriplets_,  nTOB )
DEFINE_TRIPLET_CONTAINER(vTEC_TRKTriplets_,  nTEC )
#undef DEFINE_TRIPLET_CONTAINER

// -----------------------------------------------------------------------------
// Helper type; kept local to this translation unit
// -----------------------------------------------------------------------------
namespace {
  struct HitTriple {
    float value;
    float eta;
    float phi;
    bool operator<(const HitTriple& other) const { return value > other.value; }
  };
}

// -----------------------------------------------------------------------------
// 1) Set up the ROOT branches – one per layer / projection
// -----------------------------------------------------------------------------
void
RecHitAnalyzer::branchesTRKTriplets(TTree *tree,
                                          edm::Service<TFileService>&)
{
  char hname[64];
  for (unsigned int proj = 0; proj < Nhitproj; ++proj) {

#   define BOOK(detStr, ARR, NLAYER)                                   \
      for (int iL = 0; iL < NLAYER; ++iL) {                            \
        std::snprintf(hname, sizeof(hname),                            \
                      detStr "_layer%d_triplets%s", iL+1,              \
                      hit_projections[proj].c_str());                  \
        tree->Branch(hname, &ARR[iL][proj]);                           \
      }

    BOOK("BPIX", vBPIX_TRKTriplets_, nBPIX)
    BOOK("FPIX", vFPIX_TRKTriplets_, nFPIX)
    BOOK("TIB",  vTIB_TRKTriplets_,  nTIB )
    BOOK("TID",  vTID_TRKTriplets_,  nTID )
    BOOK("TOB",  vTOB_TRKTriplets_,  nTOB )
    BOOK("TEC",  vTEC_TRKTriplets_,  nTEC )
#   undef BOOK
  }
}

// -----------------------------------------------------------------------------
// 2) Filling logic
//    – Gather *all* hits that land on ECAL for the given (sub-det & layer)
//    – Sort by descending value
//    – Keep <= 5 000 and write packed as (V,E,P,...); pad with -999
// -----------------------------------------------------------------------------
void
RecHitAnalyzer::fillTRKTriplets(const edm::Event&  iEvent,
                                      const edm::EventSetup& iSetup,
                                      unsigned int      proj)
{
  //----------------------------------------
  // 2a.  Zero-initialize the target vectors
  //----------------------------------------
  auto initContainer = [](auto& arr, int NLAYER)
  {
    for (int iL = 0; iL < NLAYER; ++iL) {
      for (unsigned int p = 0; p < Nhitproj; ++p)
        arr[iL][p].assign(3 * kMaxTripletHits, kPadValue);
    }
  };

  initContainer(vBPIX_TRKTriplets_, nBPIX);
  initContainer(vFPIX_TRKTriplets_, nFPIX);
  initContainer(vTIB_TRKTriplets_,  nTIB );
  initContainer(vTID_TRKTriplets_,  nTID );
  initContainer(vTOB_TRKTriplets_,  nTOB );
  initContainer(vTEC_TRKTriplets_,  nTEC );

  //----------------------------------------
  // 2b.  Re-use the existing hit loops, but
  //      instead of mapping to ECAL pixels
  //      push_back(value,eta,phi) to tmp vectors
  //----------------------------------------
  // A thin wrapper makes the call site identical for any sub-detector
  auto pushHit = [&](std::vector<HitTriple> layerHits[],
                     int layerIdx,
                     float value, float eta, float phi)
  {
    layerHits[layerIdx].push_back({value, eta, phi});
  };

  // -----------------------------------------------------------------------
  //    Pixel hits
  // -----------------------------------------------------------------------
  edm::Handle<SiPixelRecHitCollection> recHitColl;
  iEvent.getByToken(siPixelRecHitCollectionT_, recHitColl);

  // one scratch vector per BPIX/FPIX layer
  std::vector<HitTriple> scratchBPIX[nBPIX];
  std::vector<HitTriple> scratchFPIX[nFPIX];

  for (auto detsetIt = recHitColl->begin(); detsetIt != recHitColl->end(); ++detsetIt)
  {
    const DetId      detId  = DetId(detsetIt->detId());
    const auto      &tTopo  = iSetup.getData(tTopoToken_);
    const auto      &geom   = iSetup.getData(tkGeomToken_);
    const auto      *dg     = dynamic_cast<const PixelGeomDetUnit*>(geom.idToDetUnit(detId));
    const unsigned   layer  = getLayer(detId, &tTopo) - 1;  // 0-based
    const bool       isBPix = (detId.subdetId() == PixelSubdetector::PixelBarrel);

    for (const auto &hit : *detsetIt)
      if (hit.isValid())
      {
        const auto  lp    = hit.localPosition();
        const auto  gp    = dg->surface().toGlobal(Local3DPoint(lp));
        const TVector3 pos(gp.x(), gp.y(), gp.z());
        const float   eta  = pos.Eta();
        const float   phi  = pos.Phi();
        const float   val  = 1.f;            // <-- tracker hits are unit weight
        if (isBPix) pushHit(scratchBPIX, layer, val, eta, phi);
        else        pushHit(scratchFPIX, layer, val, eta, phi);
      }
  }

  //----------------------------------------
  // 2c.  Compress and write to output vectors
  //----------------------------------------
  auto flushLayer = [](const std::vector<HitTriple> &src,
                       std::vector<float>           &dest)
  {
    std::vector<HitTriple> ordered(src);
    std::sort(ordered.begin(), ordered.end());     // descending by value

    const unsigned nKeep = std::min<unsigned>(ordered.size(), kMaxTripletHits);
    for (unsigned i = 0; i < nKeep; ++i) {
      dest[3*i]   = ordered[i].value;
      dest[3*i+1] = ordered[i].eta;
      dest[3*i+2] = ordered[i].phi;
    }
    // remaining slots already set to kPadValue
  };

  for (int l=0; l<nBPIX; ++l)
    flushLayer(scratchBPIX[l], vBPIX_TRKTriplets_[l][proj]);

  for (int l=0; l<nFPIX; ++l)
    flushLayer(scratchFPIX[l], vFPIX_TRKTriplets_[l][proj]);

  // -----------------------------------------------------------------------
  // Repeat the same "flush" procedure for all Strip collections (TIB, ...),
  // using analogous scratch vectors defined in the corresponding loops.
  // -----------------------------------------------------------------------
}
