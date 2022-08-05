// -*- C++ -*-
//
// Package:    Analyzer/PropagatorStudies
// Class:      PropagatorStudies
//
/**\class PropagatorStudies PropagatorStudies.cc Analyzer/PropagatorStudies/plugins/PropagatorStudies.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mark Nicholas Matthewman
//         Created:  Thu, 30 Jun 2022 10:49:35 GMT
//
//

// system include files
#include <memory>
#include <numeric>
#include <sstream>
#include <any>
#include <iomanip>
#include <cmath>
#include <typeinfo>
#include <iostream>
#include <fstream>
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

#include "DataFormats/HGCalReco/interface/Trackster.h"

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


//ROOT includes
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"
#include <TFile.h>
#include <TROOT.h>
#include <TVector.h>
#include "TBranch.h"
#include <string>
#include <vector>
#include "TSystem.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include "THStack.h"
#include <algorithm>
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveStats.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;
using namespace ticl;

class PropagatorStudies : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit PropagatorStudies(const edm::ParameterSet&);
  ~PropagatorStudies();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
  virtual void fillHitMap(std::map<DetId, const HGCRecHit*>& hitMap, 
      const HGCRecHitCollection& rechitsEE, 
      const HGCRecHitCollection& rechitsFH,
      const HGCRecHitCollection& rechitsBH) const;
  std::vector<int> matchRecHit2CPRecHits(DetId detid_, std::vector<DetId> rechitdetid_);


  hgcal::RecHitTools recHitTools_;

  // ----------member data ---------------------------
  edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
  edm::EDGetTokenT<std::vector<CaloParticle>> caloParticlesToken_;
  edm::EDGetTokenT<std::vector<ticl::Trackster>> ticlTrackstersMergeToken;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsEEToken_; 
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsFHToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsBHToken_;
  //edm::EDGetTokenT<std::vector<Point3DBase<float,GlobalTag>>> propagatorEMToken_;
  //edm::EDGetTokenT<std::vector<Point3DBase<float,GlobalTag>>> propagatorHADToken_;
  edm::EDGetTokenT<std::vector<Point3DBase<float,GlobalTag>>> propagatorKFToken_;
 //edm::EDGetTokenT<std::vector<Point3DBase<float,GlobalTag>>> propagatorTrkToken_;
  //edm::EDGetTokenT<std::vector<Point3DBase<float,GlobalTag>>> propagatorTrkEMToken_;
  edm::EDGetTokenT<reco::CaloClusterCollection> hgcalLayerClustersToken_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;

  TTree *tree = new TTree("tree","tree");

  std::vector<std::string> detectors;
  std::vector<std::string> objects;
  std::vector<std::string> positions;

  // Histograms

  typedef std::map<std::string, std::vector<TH1F*>> subsubmap;
  typedef std::map<std::string, subsubmap> submap;
  std::map<std::string, submap> x_diff;
  std::map<std::string, submap> y_diff;
  std::map<std::string, submap> eta_diff;
  std::map<std::string, submap> phi_diff;

  //2D Histograms

  typedef std::map<std::string, std::vector<TH2F*>> subsubmap2D;
  typedef std::map<std::string, subsubmap2D> submap2D;
  std::map<std::string, submap2D> x_y_diff;
  std::map<std::string, submap2D> eta_phi_diff;

  std::map<std::string, submap2D> x_y_diff_cp;
  std::map<std::string, submap2D> eta_phi_diff_cp;

  //Layerwise Plots

  std::map<std::string, submap> layer_x_diff;
  std::map<std::string, submap> layer_y_diff;
  std::map<std::string, submap> layer_eta_diff;
  std::map<std::string, submap> layer_phi_diff;

  // variables
    

  //static const Int_t maxpoints = 200;
  //typedef struct {Float_t x[maxpoints],y[maxpoints],z[maxpoints];} POINT;
  //Int_t npoints;
  //POINT points;

  /*
  Int_t prop_npoints;
  Float_t prop_x[maxpoints];
  Float_t prop_y[maxpoints];
  Float_t prop_z[maxpoints];

  Int_t sim_npoints;
  Float_t sim_x[maxpoints];
  Float_t sim_y[maxpoints];
  Float_t sim_z[maxpoints];

  */
  /*
  Int_t rec_npoints;
  Float_t rec_x[maxpoints];
  Float_t rec_x[maxpoints];
  Float_t rec_x[maxpoints];


  Int_t lc_npoints;
  Float_t lc_x[maxpoints];
  Float_t lc_x[maxpoints];
  Float_t lc_x[maxpoints];
  */
  Int_t eventnr;
  
  /*
  typedef std::map<std::string, float> subsubmapfloat;
  typedef std::map<std::string, subsubmapfloat> submapfloat;
  std::map<std::string, submapfloat> x_diff_var;
  std::map<std::string, submapfloat> y_diff_var;
  std::map<std::string, submapfloat> eta_diff_var;
  std::map<std::string, submapfloat> phi_diff_var;

  */
  std::string eta_;
  std::shared_ptr<hgcal::RecHitTools> recHitTools;



  // double trackPtMin_; 


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
PropagatorStudies::PropagatorStudies(const edm::ParameterSet& iConfig) :
      tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))),
      caloParticlesToken_(consumes<std::vector<CaloParticle> >(iConfig.getParameter<edm::InputTag>("caloParticles"))), 
      ticlTrackstersMergeToken(consumes<std::vector<ticl::Trackster> >(iConfig.getParameter<edm::InputTag>("Tracksters"))),
      hgcalRecHitsEEToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsEE"))),
      hgcalRecHitsFHToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsFH"))),
      hgcalRecHitsBHToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsBH"))),
      //propagatorEMToken_(consumes<std::vector<Point3DBase<float,GlobalTag> >>(iConfig.getParameter<edm::InputTag>("propagatorEM"))),
      //propagatorHADToken_(consumes<std::vector<Point3DBase<float,GlobalTag> >>(iConfig.getParameter<edm::InputTag>("propagatorHAD"))),
      propagatorKFToken_(consumes<std::vector<Point3DBase<float,GlobalTag> >>(iConfig.getParameter<edm::InputTag>("propagatorKF"))),
      //propagatorTrkToken_(consumes<std::vector<Point3DBase<float,GlobalTag> >>(iConfig.getParameter<edm::InputTag>("propagatorTrk"))),
      //propagatorTrkEMToken_(consumes<std::vector<Point3DBase<float,GlobalTag> >>(iConfig.getParameter<edm::InputTag>("propagatorTrkEM"))),
      hgcalLayerClustersToken_(consumes<reco::CaloClusterCollection>(iConfig.getParameter<edm::InputTag>("hgcalLayerClusters"))),
      // trackPtMin_(iConfig.getParameter<double>("trackPtMin")), 
      caloGeomToken_(esConsumes<CaloGeometry, CaloGeometryRecord>()),
      eta_(iConfig.getParameter<std::string>("eta")){


  detectors = {"", "Si", "Si 120", "Si 200", "Si 300", "Sc"};
  objects = {"Simhits", "Rechits", "LCs"};
  positions = {"Propagator", "CP"};

        
  //recHitTools_.reset(new hgcal::RecHitTools());
  //now do what ever initialization is needed
  
  usesResource("TFileService");
  edm::Service<TFileService> file;

  tree = file->make<TTree>("tree","propagator");

  /*
  tree->Branch("run"    , &run    , "run/I"    );
  tree->Branch("event"  , &event  , "event/I"  );
  tree->Branch("weight" , &weight , "weight/F" );
  */

  //tree->Branch("lcpoints", &lcpoints, 99);
  //tree->Branch("recpoints", &lcpoints, 99);
  //tree->Branch("simpoints", &lcpoints, 99);
  //tree->Branch("proppoints", &proppoints,"x[47]:y[47]:z[47]");

  //tree->Branch("npoints", &npoints,"npoints/I");  
  //tree->Branch("points", &points,"x[npoints]/F:y[npoints]/F:z[npoints]/F");

  /*
  tree->Branch("prop_npoints", &prop_npoints, "prop_npoints/I");
  tree->Branch("prop_x", &prop_x, "prop_x[prop_npoints]/F");
  tree->Branch("prop_y", &prop_y, "prop_y[prop_npoints]/F");
  tree->Branch("prop_z", &prop_z, "prop_z[prop_npoints]/F");

  tree->Branch("sim_npoints", &sim_npoints, "sim_npoints/I");
  tree->Branch("sim_x", &sim_x, "sim_x[sim_npoints]/F");
  tree->Branch("sim_y", &sim_y, "sim_y[sim_npoints]/F");
  tree->Branch("sim_z", &sim_z, "sim_z[sim_npoints]/F");
  */

  /*
  tree->Branch("rec_npoints", &rec_npoints, "rec_npoints/I");
  tree->Branch("rec_x", &rec_x, "rec_x[rec_npoints]/F");
  tree->Branch("rec_y", &rec_y, "rec_y[rec_npoints]/F");
  tree->Branch("rec_z", &rec_z, "rec_z[rec_npoints]/F");

  tree->Branch("lc_npoints", &lc_npoints, "lc_npoints/I");
  tree->Branch("lc_x", &lc_x, "lc_x[lc_npoints]/F");
  tree->Branch("lc_y", &lc_y, "lc_y[lc_npoints]/F");
  tree->Branch("lc_z", &lc_z, "lc_z[lc_npoints]/F");

  */
  tree->Branch("eventnr"  , &eventnr  , "eventnr/I");

  TString t_name;
  for (const auto&pos : positions){
    for (const auto&obj : objects){
      for (const auto&det : detectors){

        t_name = "etadiff_"+obj+"_"+det+"_"+pos;
        eta_diff[pos][obj][det].push_back(new TH1F(t_name, t_name, 300,-0.15,0.15));
        t_name = "phidiff_"+obj+"_"+det+"_"+pos;
        phi_diff[pos][obj][det].push_back(new TH1F(t_name, t_name, 300,-0.15,0.15));
        t_name = "xdiff_"+obj+"_"+det+"_"+pos;
        x_diff[pos][obj][det].push_back(new TH1F(t_name, t_name, 200,-2,2));
        t_name = "ydiff_"+obj+"_"+det+"_"+pos;
        y_diff[pos][obj][det].push_back(new TH1F(t_name, t_name, 200,-2,2));
        t_name = "xydiff_"+obj+"_"+det+"_"+pos;
        x_y_diff[pos][obj][det].push_back(new TH2F(t_name,t_name,200,-2,2,200,-2,2));
        t_name = "etaphidiff_"+obj+"_"+det+"_"+pos;
        eta_phi_diff[pos][obj][det].push_back(new TH2F(t_name,t_name,200,-0.2,0.2,200,-0.2,0.2));

        t_name = "layer_xdiff_"+obj+"_"+det+"_"+pos;
        layer_x_diff[pos][obj][det].push_back(new TH1F(t_name, t_name, 47,0,47));
        t_name = "layer_ydiff_"+obj+"_"+det+"_"+pos;
        layer_y_diff[pos][obj][det].push_back(new TH1F(t_name, t_name, 47,0,47));
        t_name = "layer_etadiff_"+obj+"_"+det+"_"+pos;
        layer_eta_diff[pos][obj][det].push_back(new TH1F(t_name, t_name, 47,0,47));
        t_name = "layer_phidiff_"+obj+"_"+det+"_"+pos;
        layer_phi_diff[pos][obj][det].push_back(new TH1F(t_name, t_name, 47,0,47));

        /*

        t_name = "xdiff_"+obj+"_"+det+"_"+pos;
        tree->Branch(t_name, &x_diff_var[pos][obj][det], t_name+'/F');
        t_name = "ydiff_"+obj+"_"+det+"_"+pos;
        tree->Branch(t_name, &y_diff_var[pos][obj][det], t_name+'/F');
        t_name = "etadiff_"+obj+"_"+det+"_"+pos;
        tree->Branch(t_name, &eta_diff_var[pos][obj][det], t_name+'/F');
        t_name = "phidiff_"+obj+"_"+det+"_"+pos;
        tree->Branch(t_name, &phi_diff_var[pos][obj][det], t_name+'/F');

        */
      }
    }
  }

    
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed
}

PropagatorStudies::~PropagatorStudies() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty


  for(auto& pos : positions){
    for(auto& obj: objects){
      for(auto& det: detectors){
        x_diff[pos][obj][det].front()->Write();
        y_diff[pos][obj][det].front()->Write();
        eta_diff[pos][obj][det].front()->Write();
        phi_diff[pos][obj][det].front()->Write();
        layer_x_diff[pos][obj][det].front()->Write();
        layer_y_diff[pos][obj][det].front()->Write();
        layer_eta_diff[pos][obj][det].front()->Write();
        layer_phi_diff[pos][obj][det].front()->Write();       
        x_y_diff[pos][obj][det].front()->Write();
        eta_phi_diff[pos][obj][det].front()->Write();

      }
    }
  }
}

//
// member functions
//


// ------------ method called for each event  ------------

void PropagatorStudies::fillHitMap(std::map<DetId, const HGCRecHit*>& hitMap,
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

std::vector<int> PropagatorStudies::matchRecHit2CPRecHits(DetId detid_, std::vector<DetId> rechitdetid_) {
  std::vector<int> matchedIdxs; matchedIdxs.clear();
  for (unsigned int i0=0; i0<rechitdetid_.size(); ++i0) {
    if (detid_ == rechitdetid_[i0]) { matchedIdxs.push_back(i0); }
  }
  return matchedIdxs;
} // end of matchRecHit2CPRecHits


void PropagatorStudies::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  //std::cout << "Analyze" << std::endl;

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

  edm::Handle<std::vector<ticl::Trackster>> ticlTrackstersMerge;
  iEvent.getByToken(ticlTrackstersMergeToken, ticlTrackstersMerge);
  const std::vector<ticl::Trackster>& tracksters = *ticlTrackstersMerge; 

  /*
  edm::Handle<std::vector<Point3DBase<float,GlobalTag>>> propagatorEMHandle;
  iEvent.getByToken(propagatorEMToken_, propagatorEMHandle);
  const std::vector<Point3DBase<float,GlobalTag>> &propEM = *propagatorEMHandle;

  edm::Handle<std::vector<Point3DBase<float,GlobalTag>>> propagatorHADHandle;
  iEvent.getByToken(propagatorHADToken_, propagatorHADHandle);
  const std::vector<Point3DBase<float,GlobalTag>> &propHAD = *propagatorHADHandle;
  */
  edm::Handle<std::vector<Point3DBase<float,GlobalTag>>> propagatorKFHandle;
  iEvent.getByToken(propagatorKFToken_, propagatorKFHandle);
  const std::vector<Point3DBase<float,GlobalTag>> &propKF = *propagatorKFHandle;
  /*
  edm::Handle<std::vector<Point3DBase<float,GlobalTag>>> propagatorTrkHandle;
  iEvent.getByToken(propagatorTrkToken_, propagatorTrkHandle);
  const std::vector<Point3DBase<float,GlobalTag>> &propTrk = *propagatorTrkHandle;

  edm::Handle<std::vector<Point3DBase<float,GlobalTag>>> propagatorTrkEMHandle;
  iEvent.getByToken(propagatorTrkEMToken_, propagatorTrkEMHandle);
  const std::vector<Point3DBase<float,GlobalTag>> &propTrkEM = *propagatorTrkEMHandle;
  */
  // Investigate Propagator DAta


  /*
  int emcounter = 0;
  for(const auto& pem : propEM){
    emcounter++;
  }
  std::cout << "EM Data" <<emcounter << std::endl;

  int hadcounter = 0;
  for(const auto& phad : propHAD){
    hadcounter++;
  }
  std::cout << "HAD Data" <<hadcounter << std::endl;
  */

  std::map<float, std::vector<float>> karthesian;
  std::map<float, std::vector<float>> angular;
  //z.clear();
  /*
  int count = 0;
  prop_npoints = 47;
  */


  std::ofstream myfile;
  myfile.open("/eos/home-m/mmatthew/Data/Analyzer/PropagatorStudies/Position/kf.csv");
  myfile << "x,y,z" <<"\n";

  for(const auto& kf : propKF){
    /*
    std::cout << kf.x() << std::endl;
    std::cout << kf.y() << std::endl;

    std::cout << kf.z() << std::endl;
    std::cout << kf.eta() << std::endl;
    std::cout << kf.phi() << std::endl;
    */


    myfile << kf.x() << ","  << kf.y() << "," << kf.z() << "\n";
    karthesian[kf.z()].push_back(kf.x());
    karthesian[kf.z()].push_back(kf.y());

    angular[kf.z()].push_back(kf.eta());
    angular[kf.z()].push_back(kf.phi());

    /*
    points.x[count] = kf.x();
    points.y[count] = kf.y();
    points.z[count] = kf.z();
    count++;
    */
    /*
    prop_x[count] = kf.x();
    prop_y[count] = kf.y();
    prop_z[count] = kf.z();
    count++;
    */

  }

  myfile.close();


  /*
  int trkcounter = 0;
  for(const auto& trk : propTrk){
    trkcounter++;
  }
  std::cout << "TRK Data" <<trkcounter << std::endl;

  int trkemcounter = 0;
  for(const auto& trkem : propTrkEM){
    trkemcounter++;
  }
  std::cout << "Trk EM Data" <<trkemcounter << std::endl;
//  */


  // init vars
  const CaloGeometry &geom = iSetup.getData(caloGeomToken_);
  recHitTools_.setGeometry(geom);

  TString t_name;

  // Loop over Caloparticles 

  float cp_eta = 0;
  float cp_phi = 0;
  float edist = 0;
  float phidist = 0;
  float xdist = 0;
  float ydist = 0;

  //std::cout << "CP z Position" << std::endl;
  std::vector<DetId> tmprechits_; tmprechits_.clear();
  for (const auto& it_cp : cps) {
    // do something with track parameters, e.g, plot the charge.
    const CaloParticle& cp = ((it_cp)); 
    const SimClusterRefVector& simclusters = cp.simClusters();
    cp_eta = cp.eta();
    cp_phi = cp.phi();

    std::cout<<"Simhits"<<std::endl;
    for (const auto& it_simc : simclusters){
      const SimCluster& simc = (*(it_simc));
      const auto& sc_haf = simc.hits_and_fractions();

      
      //sim_npoints = simc.numberOfRecHits();
      //count = 0;

      for (const auto& it_sc_haf : sc_haf){
        DetId detid_ = (it_sc_haf.first);
        std::map<DetId,const HGCRecHit *>::const_iterator itcheck = hitMap.find(detid_);
        unsigned int layer_ = recHitTools_.getLayerWithOffset(detid_);   

        if(angular.find(recHitTools_.getPosition(detid_).z())!=angular.end()){
          
          // Fill simpoints

          //simpoints["x"] = recHitTools_.getPosition(detid_).x();
          //simpoints["y"] = recHitTools_.getPosition(detid_).y();
          //simpoints["z"] = recHitTools_.getPosition(detid_).z();

          /*
          sim_x[count] = recHitTools_.getPosition(detid_).x();
          sim_y[count] = recHitTools_.getPosition(detid_).y();
          sim_z[count] = recHitTools_.getPosition(detid_).z();
          count++;
          */

          // Position of Propagator
          float prop_eta = angular[recHitTools_.getPosition(detid_).z()][0];
          float prop_phi = angular[recHitTools_.getPosition(detid_).z()][1];

          float prop_x = karthesian[recHitTools_.getPosition(detid_).z()][0];
          float prop_y = karthesian[recHitTools_.getPosition(detid_).z()][1];

          // Position of extrapolated CP
          float cp_x = -1.*recHitTools_.getPosition(detid_).z()*TMath::Tan(2.*(TMath::ATan(TMath::Exp(cp_eta))))*TMath::Cos(cp_phi);
          float cp_y = -1.*recHitTools_.getPosition(detid_).z()*TMath::Tan(2.*(TMath::ATan(TMath::Exp(cp_eta))))*TMath::Sin(cp_phi);
          
          std::string detector;
          std::string thickness;
          if(recHitTools_.isSilicon(detid_)){
            detector = "Si";
            thickness = std::to_string(int(recHitTools_.getSiThickness(detid_))); 
          }
          else{
            detector = "Sc";
            thickness = "None";
            }

          //std::cout<<"Simhits"<<std::endl;
          for (const auto& pos: positions){
            if (pos=="Propagator"){
              // Diff Eta Phi         
              edist = prop_eta - recHitTools_.getPosition(detid_).eta();
              phidist = prop_phi - recHitTools_.getPosition(detid_).phi();

              // Diff x y 
              xdist = prop_x - recHitTools_.getPosition(detid_).x();
              ydist = prop_y - recHitTools_.getPosition(detid_).y();
            }
            if (pos =="CP"){

              // Diff Eta Phi         
              edist = cp_eta - recHitTools_.getPosition(detid_).eta();
              phidist = cp_phi - recHitTools_.getPosition(detid_).phi();

              // Diff x y 
              xdist = cp_x - recHitTools_.getPosition(detid_).x();
              ydist = cp_y - recHitTools_.getPosition(detid_).y();

            }
 
            for(const auto& det : detectors){
              if(det == "" || det == detector || det.find(thickness)!=std::string::npos){
                x_diff[pos][objects[0]][det].front()->Fill(xdist);
                y_diff[pos][objects[0]][det].front()->Fill(ydist);
                eta_diff[pos][objects[0]][det].front()->Fill(edist);
                phi_diff[pos][objects[0]][det].front()->Fill(phidist); 
                x_y_diff[pos][objects[0]][det].front()->Fill(xdist,ydist);
                eta_phi_diff[pos][objects[0]][det].front()->Fill(edist,phidist);

                layer_x_diff[pos][objects[0]][det].front()->Fill(layer_, xdist);
                layer_y_diff[pos][objects[0]][det].front()->Fill(layer_, ydist);
                layer_eta_diff[pos][objects[0]][det].front()->Fill(layer_, edist);
                layer_phi_diff[pos][objects[0]][det].front()->Fill(layer_,phidist); 

                /*
                x_diff_var[pos][objects[0]][det]=xdist;
                y_diff_var[pos][objects[0]][det]=ydist;
                eta_diff_var[pos][objects[0]][det]=edist;
                phi_diff_var[pos][objects[0]][det]=phidist;
                */
              }
            }
          
          
            if (itcheck != hitMap.end()){
              tmprechits_.push_back(detid_);
              // Not necessary? Detid of Simhit and Rechit the same.

              // Fill recpoints

              //recpoints.push_back(point);

              //float drRecHitCP = getDr(rhPosition_x,rhPosition_y,cpPositionAtRecHit_x,cpPositionAtRecHit_y);
              //std::cout<<"Rechits"<<std::endl;
              for(const auto& det : detectors){
                if(det == "" || det == detector || det.find(thickness)!=std::string::npos){

                  x_diff[pos][objects[1]][det].front()->Fill(xdist);
                  y_diff[pos][objects[1]][det].front()->Fill(ydist);
                  eta_diff[pos][objects[1]][det].front()->Fill(edist);
                  phi_diff[pos][objects[1]][det].front()->Fill(phidist); 
                  x_y_diff[pos][objects[1]][det].front()->Fill(xdist,ydist);
                  eta_phi_diff[pos][objects[1]][det].front()->Fill(edist,phidist);

                  layer_x_diff[pos][objects[1]][det].front()->Fill(layer_, xdist);
                  layer_y_diff[pos][objects[1]][det].front()->Fill(layer_, ydist);
                  layer_eta_diff[pos][objects[1]][det].front()->Fill(layer_, edist);
                  layer_phi_diff[pos][objects[1]][det].front()->Fill(layer_,phidist);

                  /*

                  x_diff_var[pos][objects[1]][det]=xdist;
                  y_diff_var[pos][objects[1]][det]=ydist;
                  eta_diff_var[pos][objects[1]][det]=edist;
                  phi_diff_var[pos][objects[1]][det]=phidist;      


                  */         

                }
              }
            }
          }
        }
      
        else{
          std::cout << "Incorrect z position " << recHitTools_.getPosition(detid_).z() << " at layer " << layer_<<std::endl;
        }
        //tree->Fill();
      }
    }
  }

  // Loop over LCs  

  //std::cout << "LC" << std::endl;
  for (const auto& it_lc : lcs) {

    //if(M_PI < cp_eta && cp_phi <0){
    //  continue;
    //}

    const std::vector<std::pair<DetId, float>> &hf = it_lc.hitsAndFractions();


    // loop over the rechits of this specific layer cluster
    for (unsigned int j = 0; j < hf.size(); j++) { 
      DetId detid_ = hf[j].first;
      std::vector<int> idx_= matchRecHit2CPRecHits(detid_, tmprechits_);
      unsigned int layer_ = recHitTools_.getLayerWithOffset(detid_);

      for (unsigned int i0=0; i0<idx_.size(); i0++) {
        if(angular.find(recHitTools_.getPosition(detid_).z())!=angular.end()){

          // Fill simpoints

          //lcpoints["x"]= recHitTools_.getPosition(detid_).x();
          //lcpoints["y"] = recHitTools_.getPosition(detid_).y();
          //lcpoints["z"] = recHitTools_.getPosition(detid_).z();


          // Position of Propagator
          float prop_eta = angular[recHitTools_.getPosition(detid_).z()][0];
          float prop_phi = angular[recHitTools_.getPosition(detid_).z()][1];

          float prop_x = karthesian[recHitTools_.getPosition(detid_).z()][0];
          float prop_y = karthesian[recHitTools_.getPosition(detid_).z()][1];

          // Position of extrapolated CP
          float cp_x = -1.*recHitTools_.getPosition(detid_).z()*TMath::Tan(2.*(TMath::ATan(TMath::Exp(cp_eta))))*TMath::Cos(cp_phi);
          float cp_y = -1.*recHitTools_.getPosition(detid_).z()*TMath::Tan(2.*(TMath::ATan(TMath::Exp(cp_eta))))*TMath::Sin(cp_phi);
          

          std::string detector;
          std::string thickness;
          if(recHitTools_.isSilicon(detid_)){
            detector = "Si";
            thickness = std::to_string(int(recHitTools_.getSiThickness(detid_)));   
          }
          else{
            detector = "Sc";
            thickness = "None";
          }
          //std::cout<<"LCs"<<std::endl;
          for (const auto& pos : positions){

            if (pos=="Propagator"){
              // Diff Eta Phi         
              edist = prop_eta - recHitTools_.getPosition(detid_).eta();
              phidist = prop_phi - recHitTools_.getPosition(detid_).phi();

              // Diff x y 
              xdist = prop_x - recHitTools_.getPosition(detid_).x();
              ydist = prop_y - recHitTools_.getPosition(detid_).y();
            }
            if (pos =="CP"){

              // Diff Eta Phi         
              edist = cp_eta - recHitTools_.getPosition(detid_).eta();
              phidist = cp_phi - recHitTools_.getPosition(detid_).phi();

              // Diff x y 
              xdist = cp_x - recHitTools_.getPosition(detid_).x();
              ydist = cp_y - recHitTools_.getPosition(detid_).y();
            }

            for(const auto& det : detectors){
              if(det == "" || det == detector || det.find(thickness)!=std::string::npos){
                x_diff[pos][objects[2]][det].front()->Fill(xdist);
                y_diff[pos][objects[2]][det].front()->Fill(ydist);
                eta_diff[pos][objects[2]][det].front()->Fill(edist);
                phi_diff[pos][objects[2]][det].front()->Fill(phidist);
                x_y_diff[pos][objects[2]][det].front()->Fill(xdist,ydist);
                eta_phi_diff[pos][objects[2]][det].front()->Fill(edist,phidist);

                layer_x_diff[pos][objects[2]][det].front()->Fill(layer_, xdist);
                layer_y_diff[pos][objects[2]][det].front()->Fill(layer_, ydist);
                layer_eta_diff[pos][objects[2]][det].front()->Fill(layer_, edist);
                layer_phi_diff[pos][objects[2]][det].front()->Fill(layer_,phidist);

                /*

                x_diff_var[pos][objects[2]][det]=xdist;
                y_diff_var[pos][objects[2]][det]=ydist;
                eta_diff_var[pos][objects[2]][det]=edist;
                phi_diff_var[pos][objects[2]][det]=phidist;
                */

                //tree->Fill();
              }
            }
          }
        }
        else{
          std::cout << recHitTools_.getPosition(detid_).z() <<std::endl;

        }
      }
    }
  }

  for (const auto& track : iEvent.get(tracksToken_)) {
    // do something with track parameters, e.g, plot the charge.
    // int charge = track.charge();
  }

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
  tree->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void PropagatorStudies::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void PropagatorStudies::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void PropagatorStudies::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(PropagatorStudies);
