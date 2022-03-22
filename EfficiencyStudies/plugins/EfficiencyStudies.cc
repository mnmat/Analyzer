// -*- C++ -*-
//
// Package:    Demo/EfficiencyStudies
// Class:      EfficiencyStudies
//
/**\class EfficiencyStudies EfficiencyStudies.cc Demo/EfficiencyStudies/plugins/EfficiencyStudies.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Shengquan Tuo
//         Created:  Thu, 23 Sep 2021 21:03:52 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"

#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticleFwd.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"

#include "SimDataFormats/CaloAnalysis/interface/SimClusterFwd.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"

#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"

#include "DataFormats/ParticleFlowReco/interface/HGCalMultiCluster.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "Geometry/CaloTopology/interface/HGCalTopology.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"



//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;

class EfficiencyStudies : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit EfficiencyStudies(const edm::ParameterSet&);
  ~EfficiencyStudies();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
  virtual void fillHitMap(std::map<DetId, const HGCRecHit*>& hitMap, 
			  const HGCRecHitCollection& rechitsEE, 
			  const HGCRecHitCollection& rechitsFH,
			  const HGCRecHitCollection& rechitsBH) const;

  hgcal::RecHitTools recHitTools_;

  // ----------member data ---------------------------
  edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
  edm::EDGetTokenT<std::vector<CaloParticle>> caloParticlesToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsEEToken_; 
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsFHToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsBHToken_;
  edm::EDGetTokenT<reco::CaloClusterCollection> hgcalLayerClustersToken_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;

  double trackPtMin_;

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
EfficiencyStudies::EfficiencyStudies(const edm::ParameterSet& iConfig) :
      tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))),
      caloParticlesToken_(consumes<std::vector<CaloParticle> >(iConfig.getParameter<edm::InputTag>("caloParticles"))), 
      hgcalRecHitsEEToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsEE"))),
      hgcalRecHitsFHToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsFH"))),
      hgcalRecHitsBHToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsBH"))),
      hgcalLayerClustersToken_(consumes<reco::CaloClusterCollection>(iConfig.getParameter<edm::InputTag>("hgcalLayerClusters"))),
      caloGeomToken_(esConsumes<CaloGeometry, CaloGeometryRecord>()),
      trackPtMin_(iConfig.getParameter<double>("trackPtMin")) {


#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif

  //recHitTools.reset(new hgcal::RecHitTools());
  //now do what ever initialization is needed
}

EfficiencyStudies::~EfficiencyStudies() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------

void EfficiencyStudies::fillHitMap(std::map<DetId, const HGCRecHit*>& hitMap,
                                const HGCRecHitCollection& rechitsEE,
                                const HGCRecHitCollection& rechitsFH,
                                const HGCRecHitCollection& rechitsBH) const {
  hitMap.clear();
  for (const auto& hit : rechitsEE) {
    hitMap.emplace(hit.detid(), &hit);
  }

  for (const auto& hit : rechitsFH) {
    hitMap.emplace(hit.detid(), &hit);
  }

  for (const auto& hit : rechitsBH) {
    hitMap.emplace(hit.detid(), &hit);
  }
} // end of EfficiencyStudies::fillHitMap



void EfficiencyStudies::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace std;

  edm::Handle<std::vector<CaloParticle>> CaloParticles;
  iEvent.getByToken(caloParticlesToken_, CaloParticles);
  const CaloParticleCollection& cps = *CaloParticles;

  
  edm::Handle<HGCRecHitCollection> recHitHandleEE;  
  iEvent.getByToken(hgcalRecHitsEEToken_, recHitHandleEE);

  edm::Handle<HGCRecHitCollection> recHitHandleFH;
  iEvent.getByToken(hgcalRecHitsFHToken_, recHitHandleFH);

  edm::Handle<HGCRecHitCollection> recHitHandleBH;
  iEvent.getByToken(hgcalRecHitsBHToken_, recHitHandleBH);

  edm::Handle<reco::CaloClusterCollection> layerClusterHandle;
  iEvent.getByToken(hgcalLayerClustersToken_, layerClusterHandle);
  const reco::CaloClusterCollection &lcs = *layerClusterHandle;

  std::map<DetId, const HGCRecHit*> hitMap;
  fillHitMap(hitMap, *recHitHandleEE, *recHitHandleFH, *recHitHandleBH);

  // init vars
  const CaloGeometry &geom = iSetup.getData(caloGeomToken_);
  //edm::ESHandle<CaloGeometry> geom = iSetup.getHandle(caloGeomToken_);
  //edm::ESHandle<CaloGeometry> geom; 
  //edm::ESGetToken<CaloGeometry,CaloGeometryRecord> geomToken_; 
  //iSetup.get<CaloGeometryRecord>().get(geom); 
  recHitTools_.setGeometry(geom);


  int nTrack = 0;
  for (const auto& track : iEvent.get(tracksToken_)) {
    // do something with track parameters, e.g, plot the charge.
    // int charge = track.charge();
    if(track.pt() < trackPtMin_) continue;
    nTrack++;
  }

  int nCalo = 0;
  int nSimhits = 0;
  int nRechits = 0;
  int nLChits = 0;
  for (const auto& it_cp : cps) {
    // do something with track parameters, e.g, plot the charge.
    // int charge = track.charge();
    nCalo++;
    const CaloParticle& cp = ((it_cp)); 

    const SimClusterRefVector& simclusters = cp.simClusters();
    for (const auto& it_simc : simclusters){

      const SimCluster& simc = (*(it_simc));
      const auto& sc_haf = simc.hits_and_fractions();

      for (const auto& it_sc_haf : sc_haf){
        nSimhits++;
        DetId detid_ = (it_sc_haf.first);
        std::map<DetId,const HGCRecHit *>::const_iterator itcheck = hitMap.find(detid_);
	      unsigned int layer_ = recHitTools_.getLayerWithOffset(detid_); 
        std::cout << "layer = " << layer_ << "\n";                                                                                      

        if (itcheck != hitMap.end()){
          nRechits++;
          }
        //cout<<"DetId = "<<detid_<<endl;         
        }
     }
  }
  
  
  for (const auto& it_lc : lcs) {
    const std::vector<std::pair<DetId, float>> &hf = it_lc.hitsAndFractions();
    // loop over the rechits of this specific layer cluster
    for (unsigned int j = 0; j < hf.size(); j++) { 
      nLChits++;
      }
    }
  
  
  cout<<"nTrack = "<<nTrack<<endl;
  cout<<"nCalo = "<<nCalo<<endl;
  cout<<"nSim = "<<nSimhits<<endl;
  cout<<"nRec = "<<nRechits<<endl;
  cout << "Num of LCs = "<< lcs.size() <<endl ;
  cout<<"nLChits = "<<nLChits<<endl;



#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void EfficiencyStudies::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void EfficiencyStudies::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void EfficiencyStudies::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(EfficiencyStudies);
