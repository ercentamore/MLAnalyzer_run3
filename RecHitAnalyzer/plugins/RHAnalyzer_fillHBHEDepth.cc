#include "MLAnalyzer_run3/RecHitAnalyzer/interface/RecHitAnalyzer.h"

static constexpr int HBHE_MAX_DEPTH = 7;
static constexpr int N_IETA_BINS = 2 * (HBHE_IETA_MAX_HE-1);
static constexpr int N_PIXELS_2D = HBHE_IPHI_NUM * N_IETA_BINS;
static constexpr int N_PIXELS_3D = HBHE_MAX_DEPTH * N_PIXELS_2D;

TH2F *hEvt_HBHE_energy;
TProfile2D *hHBHE_energy_EB;
TProfile2D *hHBHE_energy;

std::vector<float> vHBHE_energy_EB_;
std::vector<float> vHBHE_energy_;
std::vector<float> vHBHE_energy_depth_;

void RecHitAnalyzer::branchesHBHE ( TTree *tree, edm::Service<TFileService> &fs ) {
  tree->Branch( "HBHE_energy_EB", &vHBHE_energy_EB_ );
  tree->Branch( "HBHE_energy"   , &vHBHE_energy_    );

  tree->Branch( "HBHE_energy_depth", &vHBHE_energy_depth_ );

  hEvt_HBHE_energy = new TH2F("evt_HBHE_energy",
                              "E(i#phi,i#eta);i#phi;i#eta",
                              HBHE_IPHI_NUM,           HBHE_IPHI_MIN-1, HBHE_IPHI_MAX,
                              N_IETA_BINS,            -(HBHE_IETA_MAX_HE-1), HBHE_IETA_MAX_HE-1 );

  hHBHE_energy = fs->make<TProfile2D>( "HBHE_energy",
                                       "E(i#phi,i#eta);i#phi;i#eta",
                                       HBHE_IPHI_NUM, HBHE_IPHI_MIN-1, HBHE_IPHI_MAX,
                                       2*HBHE_IETA_MAX_HE, -HBHE_IETA_MAX_HE, HBHE_IETA_MAX_HE );

  hHBHE_energy_EB = fs->make<TProfile2D>( "HBHE_energy_EB",
                                          "E(i#phi,i#eta);i#phi;i#eta",
                                          HBHE_IPHI_NUM, HBHE_IPHI_MIN-1, HBHE_IPHI_MAX,
                                          2*HBHE_IETA_MAX_EB, -HBHE_IETA_MAX_EB, HBHE_IETA_MAX_EB );
} // branchesHBHE



void RecHitAnalyzer::fillHBHE ( const edm::Event &iEvent, const edm::EventSetup &iSetup ) {
  vHBHE_energy_EB_.assign ( 2 * HBHE_IPHI_NUM * HBHE_IETA_MAX_EB         , 0.f );
  vHBHE_energy_   .assign ( N_PIXELS_2D                                    , 0.f );
  vHBHE_energy_depth_.assign( N_PIXELS_3D                                  , 0.f );
  hEvt_HBHE_energy->Reset();

  edm::Handle<HBHERecHitCollection> HBHERecHitsH_;
  iEvent.getByToken( HBHERecHitCollectionT_, HBHERecHitsH_ );

  auto addToDepthImage = [&](int depthIdx, int ietaVal, int iphiBin, float e)
  {
    if ( depthIdx < 0 || depthIdx >= HBHE_MAX_DEPTH ) return;

    iphiBin = ( (iphiBin % HBHE_IPHI_NUM) + HBHE_IPHI_NUM ) % HBHE_IPHI_NUM;

    int ietaBin = ietaVal + (HBHE_IETA_MAX_HE-1);
    if ( ietaBin < 0 || ietaBin >= N_IETA_BINS ) return;

    int idx = depthIdx * N_PIXELS_2D + ietaBin * HBHE_IPHI_NUM + iphiBin;
    vHBHE_energy_depth_[ idx ] += e;
  };

  for ( const auto & rHit : *HBHERecHitsH_ ) {

    float energy = rHit.energy();
    if ( energy <= zs ) continue;

    HcalDetId hId( rHit.id() );

    int iphi  = hId.iphi() + 2;
    iphi      = iphi > HBHE_IPHI_MAX ? iphi-HBHE_IPHI_MAX : iphi;
    --iphi;

    int ietaAbs = hId.ietaAbs() == HBHE_IETA_MAX_HE ? HBHE_IETA_MAX_HE-1
                                                    : hId.ietaAbs();
    int ieta    = hId.zside() > 0 ?  ietaAbs-1 : -ietaAbs;

    int depthIdx = hId.depth()-1;

    if ( hId.ietaAbs() > HBHE_IETA_MAX_FINE ) {

      float eHalf = 0.5f * energy;

      hHBHE_energy->Fill( iphi   , ieta , eHalf );
      hHBHE_energy->Fill( iphi+1 , ieta , eHalf );

      hEvt_HBHE_energy->Fill( iphi   , ieta , eHalf );
      hEvt_HBHE_energy->Fill( iphi+1 , ieta , eHalf );

      addToDepthImage( depthIdx, ieta,  iphi   , eHalf );
      addToDepthImage( depthIdx, ieta,  iphi+1 , eHalf );
    }

    else {

      hHBHE_energy->Fill ( iphi , ieta , energy );
      hEvt_HBHE_energy->Fill( iphi , ieta , energy );

      addToDepthImage( depthIdx, ieta, iphi, energy );
    }

    if ( hId.ietaAbs() > HBHE_IETA_MAX_EB ) continue;

    hHBHE_energy_EB->Fill( iphi , ieta , energy );

    int idxEB = ( ieta + HBHE_IETA_MAX_EB ) * HBHE_IPHI_NUM + iphi;
    vHBHE_energy_EB_[ idxEB ] += energy;
  }

  for ( int iy = 1; iy <= hEvt_HBHE_energy->GetNbinsY(); ++iy )
    for ( int ix = 1; ix <= hEvt_HBHE_energy->GetNbinsX(); ++ix ) {

      float e = hEvt_HBHE_energy->GetBinContent( ix, iy );
      if ( e <= zs ) continue;

      int idx2D = (iy-1) * HBHE_IPHI_NUM + (ix-1);
      vHBHE_energy_[ idx2D ] = e;
    }
} // fillHBHE()