// -*- C++ -*-
//
// Package:    Analyzer/UpdatorStudies
// Class:      UpdatorStudies
//
/**\class UpdatorStudies UpdatorStudies.cc Analyzer/UpdatorStudies/plugins/UpdatorStudies.cc

 Description: [one line class summary]

 Implementation:f
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
#include "TProfile.h"

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

class UpdatorStudies : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit UpdatorStudies(const edm::ParameterSet&);
  ~UpdatorStudies();

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
  bool checkHex(float x, float y, DetId detid_);
  bool checkScint(float eta, float phi, DetId detid_);


  hgcal::RecHitTools recHitTools_;

  // ----------member data ---------------------------
  edm::EDGetTokenT<std::vector<CaloParticle>> caloParticlesToken_;
  edm::EDGetTokenT<std::vector<ticl::Trackster>> ticlTrackstersMergeToken;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsEEToken_; 
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsFHToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsBHToken_;
  //edm::EDGetTokenT<std::vector<Point3DBase<float,GlobalTag>>> propagatorEMToken_;
  //edm::EDGetTokenT<std::vector<Point3DBase<float,GlobalTag>>> propagatorHADToken_;
  edm::EDGetTokenT<std::vector<Point3DBase<float,GlobalTag>>> propagatorKFToken_;
  edm::EDGetTokenT<std::vector<float>> xxKFToken_;
  edm::EDGetTokenT<std::vector<float>> xyKFToken_;
  edm::EDGetTokenT<std::vector<float>> yyKFToken_;
  edm::EDGetTokenT<std::vector<Point3DBase<float,GlobalTag>>> propagatorToken_;
  edm::EDGetTokenT<std::vector<float>> xxPropToken_;
  edm::EDGetTokenT<std::vector<float>> xyPropToken_;
  edm::EDGetTokenT<std::vector<float>> yyPropToken_;
  edm::EDGetTokenT<float> abs_failToken_;
  edm::EDGetTokenT<std::vector<float>> chargeToken_;
 //edm::EDGetTokenT<std::vector<Point3DBase<float,GlobalTag>>> propagatorTrkToken_;
  //edm::EDGetTokenT<std::vector<Point3DBase<float,GlobalTag>>> propagatorTrkEMToken_;
  edm::EDGetTokenT<reco::CaloClusterCollection> hgcalLayerClustersToken_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;

  TTree *tree = new TTree("tree","tree");

  std::vector<std::string> detectors, objects, positions, hittypes;

  // Histograms

  typedef std::map<std::string, std::vector<TH1F*>> subsubmap;
  typedef std::map<std::string, subsubmap> submap;
  std::map<std::string, submap> x_diff;
  std::map<std::string, submap> y_diff;
  std::map<std::string, submap> x_pull;
  std::map<std::string, submap> y_pull;
  std::map<std::string, submap> eta_diff;
  std::map<std::string, submap> phi_diff;
  std::map<std::string, submap> r_diff;

  int eventidx = 0;
  //2D Histograms

  typedef std::map<std::string, std::vector<TH2F*>> subsubmap2D;
  typedef std::map<std::string, subsubmap2D> submap2D;
  std::map<std::string, submap2D> x_y_diff;
  std::map<std::string, submap2D> eta_phi_diff;
  std::map<std::string, submap2D> layer_xydiff;

  //Layerwise Plots

  std::map<std::string, submap> layer_x_diff;
  std::map<std::string, submap> layer_y_diff;
  std::map<std::string, submap> layer_eta_diff;
  std::map<std::string, submap> layer_phi_diff;
  std::map<std::string, submap> layer_r_diff;

  std::map<std::string, submap> layer_abs_x_diff;
  std::map<std::string, submap> layer_abs_y_diff;
  std::map<std::string, submap> layer_abs_eta_diff;
  std::map<std::string, submap> layer_abs_phi_diff;

  std::map<std::string, submap> layer_eff_kf_prop;
  std::map<std::string, submap> layer_eff_kf_hit;
  std::map<std::string, submap> layer_eff_prop_hit;

  std::map<std::string, submap> layer_eff;

  submap layer_hits;
  submap layer_energies;
  submap hits;
  submap energies;
  submap single_layer_energies;
  submap single_layer_hits;

  //2D Profiles

  typedef std::map<std::string, std::vector<TProfile*>> subsubmapProfile;
  typedef std::map<std::string, subsubmapProfile> submapProfile;

  // 2D Histos

  std::map<std::string, submapProfile> layer_profile_abs_x_diff;
  std::map<std::string, submapProfile> layer_profile_abs_y_diff;
  std::map<std::string, submapProfile> layer_profile_abs_eta_diff;
  std::map<std::string, submapProfile> layer_profile_abs_phi_diff;

  // variables
  
  Int_t eventnr;
  std::string eta_;
  std::string energy_;
  std::string outdir_;
  std::shared_ptr<hgcal::RecHitTools> recHitTools;

  int pcharge;
  int fail;
  float hit_energy, kf_energy, prop_energy;

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
UpdatorStudies::UpdatorStudies(const edm::ParameterSet& iConfig) :
      caloParticlesToken_(consumes<std::vector<CaloParticle> >(iConfig.getParameter<edm::InputTag>("caloParticles"))), 
      ticlTrackstersMergeToken(consumes<std::vector<ticl::Trackster> >(iConfig.getParameter<edm::InputTag>("Tracksters"))),
      hgcalRecHitsEEToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsEE"))),
      hgcalRecHitsFHToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsFH"))),
      hgcalRecHitsBHToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsBH"))),
      //propagatorEMToken_(consumes<std::vector<Point3DBase<float,GlobalTag> >>(iConfig.getParameter<edm::InputTag>("propagatorEM"))),
      //propagatorHADToken_(consumes<std::vector<Point3DBase<float,GlobalTag> >>(iConfig.getParameter<edm::InputTag>("propagatorHAD"))),
      propagatorKFToken_(consumes<std::vector<Point3DBase<float,GlobalTag> >>(iConfig.getParameter<edm::InputTag>("propagatorKF"))),
      xxKFToken_(consumes<std::vector<float>>(iConfig.getParameter<edm::InputTag>("xxKF"))),
      xyKFToken_(consumes<std::vector<float>>(iConfig.getParameter<edm::InputTag>("xyKF"))),
      yyKFToken_(consumes<std::vector<float>>(iConfig.getParameter<edm::InputTag>("yyKF"))),
      propagatorToken_(consumes<std::vector<Point3DBase<float,GlobalTag> >>(iConfig.getParameter<edm::InputTag>("propagator"))),
      xxPropToken_(consumes<std::vector<float>>(iConfig.getParameter<edm::InputTag>("xxProp"))),
      xyPropToken_(consumes<std::vector<float>>(iConfig.getParameter<edm::InputTag>("xyProp"))),
      yyPropToken_(consumes<std::vector<float>>(iConfig.getParameter<edm::InputTag>("yyProp"))),
      abs_failToken_(consumes<float>(iConfig.getParameter<edm::InputTag>("abs_fail"))),
      chargeToken_(consumes<std::vector<float>>(iConfig.getParameter<edm::InputTag>("charge"))),
      //propagatorTrkToken_(consumes<std::vector<Point3DBase<float,GlobalTag> >>(iConfig.getParameter<edm::InputTag>("propagatorTrk"))),
      //propagatorTrkEMToken_(consumes<std::vector<Point3DBase<float,GlobalTag> >>(iConfig.getParameter<edm::InputTag>("propagatorTrkEM"))),
      hgcalLayerClustersToken_(consumes<reco::CaloClusterCollection>(iConfig.getParameter<edm::InputTag>("hgcalLayerClusters"))),
      caloGeomToken_(esConsumes<CaloGeometry, CaloGeometryRecord>()),
      eta_(iConfig.getParameter<std::string>("eta")),
      energy_(iConfig.getParameter<std::string>("energy")),
      outdir_(iConfig.getParameter<std::string>("outdir")){


  
  std::ofstream kffile(outdir_ + "/Positions/kf_en"+energy_+"_eta_"+eta_+".csv", std::ios::out);
  kffile << "Event"<< ","<< "x"<< ","<<"y"<< ","<<"z"<<std::endl;
  kffile.close();
  std::ofstream propfile(outdir_ + "/Positions/prop_en"+energy_+"_eta_"+eta_+".csv", std::ios::out);
  propfile << "Event"<< ","<< "x"<< ","<<"y"<< ","<<"z"<<std::endl;
  propfile.close();
  std::ofstream simfile(outdir_ + "/Positions/sim_en"+energy_+"_eta_"+eta_+".csv", std::ios::out);
  simfile << "Event"<< ","<< "x"<< ","<<"y"<< ","<<"z"<<std::endl;
  simfile.close();
  std::ofstream recfile(outdir_ + "/Positions/rec_en"+energy_+"_eta_"+eta_+".csv", std::ios::out);
  recfile << "Event"<< ","<< "x"<< ","<<"y"<< ","<<"z"<<std::endl;
  recfile.close();
  std::ofstream efffile(outdir_ + "/Efficiency/eff_en"+energy_+"_eta_"+eta_+".csv", std::ios::out);
  efffile << "Event,eff_kf,eff_prop,simhits,rechits" <<std::endl;
  efffile.close();


  detectors = {"", "Si", "Si 120", "Si 200", "Si 300", "Sc"};
  hittypes = {"Simhits","Rechits","Propagator","KF"};
  objects = {"Simhits", "Rechits"};
  positions = {"Propagator", "KF"};
        
  //recHitTools_.reset(new hgcal::RecHitTools());
  //now do what ever initialization is needed
  
  usesResource("TFileService");
  edm::Service<TFileService> file;

  tree = file->make<TTree>("tree","tree");

  /*
  tree->Branch("run"    , &run    , "run/I"    );
  tree->Branch("event"  , &event  , "event/I"  );
  tree->Branch("weight" , &weight , "weight/F" );
  */

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

  */

  tree->Branch("eventnr"  , &eventnr  , "eventnr/I");
  tree->Branch("pcharge", &pcharge);
  tree->Branch("fail", &fail);

  TString t_name;


  for (const auto&det : detectors){
    for (const auto&h : hittypes){
      t_name = "layer_hits_"+h+"_"+det;
      layer_hits[h][det].push_back(new TH1F(t_name, t_name, 47,0,47));
      t_name = "hits_"+h+"_"+det;
      hits[h][det].push_back(new TH1F(t_name, t_name, 100,0,100));
      t_name = "layer_energies_"+h+"_"+det;
      layer_energies[h][det].push_back(new TH1F(t_name, t_name, 47,0,47));
      t_name = "energies_"+h+"_"+det;
      energies[h][det].push_back(new TH1F(t_name, t_name, 100,0,1));
      for (int layer=48; layer>=0; layer--){
        std::stringstream nlayer;
        nlayer << layer;
        t_name = nlayer.str() + "_layer_energies_"+h+"_"+det;    
        single_layer_energies[h][det].push_back(new TH1F(t_name,t_name,100,0,1));
        t_name = nlayer.str() + "_layer_hits_"+h+"_"+det;    
        single_layer_hits[h][det].push_back(new TH1F(t_name,t_name,10,0,10));
      }
    }

    for (const auto&pos : positions){
      for (const auto&obj : objects){
        t_name = "etadiff_"+obj+"_"+det+"_"+pos;
        eta_diff[pos][obj][det].push_back(new TH1F(t_name, t_name, 300,-0.15,0.15));
        t_name = "phidiff_"+obj+"_"+det+"_"+pos;
        phi_diff[pos][obj][det].push_back(new TH1F(t_name, t_name, 300,-0.15,0.15));
        t_name = "xdiff_"+obj+"_"+det+"_"+pos;
        x_diff[pos][obj][det].push_back(new TH1F(t_name, t_name, 600,-3,3));
        t_name = "ydiff_"+obj+"_"+det+"_"+pos;
        y_diff[pos][obj][det].push_back(new TH1F(t_name, t_name, 600,-3,3));
        t_name = "xydiff_"+obj+"_"+det+"_"+pos;
        x_y_diff[pos][obj][det].push_back(new TH2F(t_name,t_name,600,-3,3,120,-3,3));
        t_name = "etaphidiff_"+obj+"_"+det+"_"+pos;
        eta_phi_diff[pos][obj][det].push_back(new TH2F(t_name,t_name,600,-0.2,0.2,120,-0.2,0.2));
        t_name = "rdiff_"+obj+"_"+det+"_"+pos;
        r_diff[pos][obj][det].push_back(new TH1F(t_name, t_name, 600,-3,3));
        t_name = "xpull_"+obj+"_"+det+"_"+pos;
        x_pull[pos][obj][det].push_back(new TH1F(t_name, t_name, 1400,-7,7));
        t_name = "ypull_"+obj+"_"+det+"_"+pos;
        y_pull[pos][obj][det].push_back(new TH1F(t_name, t_name, 1400,-7,7));

        t_name = "layer_xdiff_"+obj+"_"+det+"_"+pos;
        layer_x_diff[pos][obj][det].push_back(new TH1F(t_name, t_name, 47,0,47));
        t_name = "layer_ydiff_"+obj+"_"+det+"_"+pos;
        layer_y_diff[pos][obj][det].push_back(new TH1F(t_name, t_name, 47,0,47));
        t_name = "layer_etadiff_"+obj+"_"+det+"_"+pos;
        layer_eta_diff[pos][obj][det].push_back(new TH1F(t_name, t_name, 47,0,47));
        t_name = "layer_phidiff_"+obj+"_"+det+"_"+pos;
        layer_phi_diff[pos][obj][det].push_back(new TH1F(t_name, t_name, 47,0,47));        
        t_name = "layer_rdiff_"+obj+"_"+det+"_"+pos;
        layer_r_diff[pos][obj][det].push_back(new TH1F(t_name, t_name, 47,0,47));

        t_name = "layer_abs_xdiff_"+obj+"_"+det+"_"+pos;
        layer_abs_x_diff[pos][obj][det].push_back(new TH1F(t_name, t_name, 47,0,47));
        t_name = "layer_abs_ydiff_"+obj+"_"+det+"_"+pos;
        layer_abs_y_diff[pos][obj][det].push_back(new TH1F(t_name, t_name, 47,0,47));
        t_name = "layer_abs_etadiff_"+obj+"_"+det+"_"+pos;
        layer_abs_eta_diff[pos][obj][det].push_back(new TH1F(t_name, t_name, 47,0,47));
        t_name = "layer_abs_phidiff_"+obj+"_"+det+"_"+pos;
        layer_abs_phi_diff[pos][obj][det].push_back(new TH1F(t_name, t_name, 47,0,47));

        t_name = "layer_profile_abs_xdiff_"+obj+"_"+det+"_"+pos;
        layer_profile_abs_x_diff[pos][obj][det].push_back(new TProfile(t_name, t_name, 47,0,47,0,5));
        t_name = "layer_profile_abs_ydiff_"+obj+"_"+det+"_"+pos;
        layer_profile_abs_y_diff[pos][obj][det].push_back(new TProfile(t_name, t_name, 47,0,47,0,5));
        t_name = "layer_profile_abs_etadiff_"+obj+"_"+det+"_"+pos;
        layer_profile_abs_eta_diff[pos][obj][det].push_back(new TProfile(t_name, t_name, 47,0,47,0,5));
        t_name = "layer_profile_abs_phidiff_"+obj+"_"+det+"_"+pos;
        layer_profile_abs_phi_diff[pos][obj][det].push_back(new TProfile(t_name, t_name, 47,0,47,0,5));


        t_name = "layer_eff_"+obj+"_"+det+"_"+pos;
        layer_eff[pos][obj][det].push_back(new TH1F(t_name, t_name, 47,0,47));

        std::map<std::string, submap> layer_eff_kf_prop;
        std::map<std::string, submap> layer_eff_kf_hit;
        std::map<std::string, submap> layer_eff_prop_hit;

        for(int i=0; i<48; i++){
          std::stringstream nlayer;
          nlayer << i;
          t_name = "layer_xydiff_" + obj + "_" + det + "_" + pos +"_Layer_" + nlayer.str();    
          layer_xydiff[pos][obj][det].push_back(new TH2F(t_name,t_name,600,-3,3,600,-3,3));
        }
      }
    }
  }

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed
}

UpdatorStudies::~UpdatorStudies() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty

  for(auto& det: detectors){
    for(auto&h: hittypes){
      layer_hits[h][det].front()->Write();
      hits[h][det].front()->Write();
      layer_energies[h][det].front()->Write();
      energies[h][det].front()->Write();
      for (int layer=0;layer<48; layer++){
        single_layer_energies[h][det][layer]->Write();
        single_layer_hits[h][det][layer]->Write();
      }
    }
    for(auto& pos : positions){
      for(auto& obj: objects){
        x_diff[pos][obj][det].front()->Write();
        y_diff[pos][obj][det].front()->Write();
        eta_diff[pos][obj][det].front()->Write();
        phi_diff[pos][obj][det].front()->Write();
        r_diff[pos][obj][det].front()->Write();
        x_pull[pos][obj][det].front()->Write();
        y_pull[pos][obj][det].front()->Write();
        layer_x_diff[pos][obj][det].front()->Write();
        layer_y_diff[pos][obj][det].front()->Write();
        layer_eta_diff[pos][obj][det].front()->Write();
        layer_phi_diff[pos][obj][det].front()->Write();
        layer_abs_x_diff[pos][obj][det].front()->Write();
        layer_abs_y_diff[pos][obj][det].front()->Write();
        layer_abs_eta_diff[pos][obj][det].front()->Write();
        layer_abs_phi_diff[pos][obj][det].front()->Write();   
        layer_profile_abs_x_diff[pos][obj][det].front()->Write();
        layer_profile_abs_y_diff[pos][obj][det].front()->Write();
        layer_profile_abs_eta_diff[pos][obj][det].front()->Write();
        layer_profile_abs_phi_diff[pos][obj][det].front()->Write(); 
        x_y_diff[pos][obj][det].front()->Write();
        eta_phi_diff[pos][obj][det].front()->Write();
        layer_eff[pos][obj][det].front()->Write();

        for(int i=0; i<48; i++){
          layer_xydiff[pos][obj][det][i]->Write();
        }
      }
    }
  }
}

//
// member functions
//

// ------------ method called for each event  ------------

void UpdatorStudies::fillHitMap(std::map<DetId, const HGCRecHit*>& hitMap,
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


std::vector<int> UpdatorStudies::matchRecHit2CPRecHits(DetId detid_, std::vector<DetId> rechitdetid_) {
  std::vector<int> matchedIdxs; matchedIdxs.clear();
  for (unsigned int i0=0; i0<rechitdetid_.size(); ++i0) {
    if (detid_ == rechitdetid_[i0]) { matchedIdxs.push_back(i0); }
  }
  return matchedIdxs;
} // end of matchRecHit2CPRecHits


void UpdatorStudies::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  eventidx +=1;
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

  std::map<DetId, const HGCRecHit*> hitMap;
  fillHitMap(hitMap, *recHitHandleEE, *recHitHandleFH, *recHitHandleBH);

  edm::Handle<std::vector<ticl::Trackster>> ticlTrackstersMerge;
  iEvent.getByToken(ticlTrackstersMergeToken, ticlTrackstersMerge);
  const std::vector<ticl::Trackster>& tracksters = *ticlTrackstersMerge; 

  edm::Handle<std::vector<Point3DBase<float,GlobalTag>>> propagatorKFHandle;
  iEvent.getByToken(propagatorKFToken_, propagatorKFHandle);
  const std::vector<Point3DBase<float,GlobalTag>> &kf = *propagatorKFHandle;

  edm::Handle<std::vector<Point3DBase<float,GlobalTag>>> propagatorHandle;
  iEvent.getByToken(propagatorToken_, propagatorHandle);
  const std::vector<Point3DBase<float,GlobalTag>> &prop = *propagatorHandle;

  edm::Handle<std::vector<float>> xxPropHandle;
  iEvent.getByToken(xxPropToken_, xxPropHandle);
  const std::vector<float> &xxprop = *xxPropHandle;

  edm::Handle<std::vector<float>> xyPropHandle;
  iEvent.getByToken(xyPropToken_, xyPropHandle);
  const std::vector<float> &xyprop = *xyPropHandle;

  edm::Handle<std::vector<float>> yyPropHandle;
  iEvent.getByToken(yyPropToken_, yyPropHandle);
  const std::vector<float> &yyprop = *yyPropHandle;

  edm::Handle<std::vector<float>> xxKFHandle;
  iEvent.getByToken(xxKFToken_, xxKFHandle);
  const std::vector<float> &xxkf = *xxKFHandle;

  edm::Handle<std::vector<float>> xyKFHandle;
  iEvent.getByToken(xyKFToken_, xyKFHandle);
  const std::vector<float> &xykf = *xyKFHandle;

  edm::Handle<std::vector<float>> yyKFHandle;
  iEvent.getByToken(yyKFToken_, yyKFHandle);
  const std::vector<float> &yykf = *yyKFHandle;

  edm::Handle<float> abs_failHandle_;
  iEvent.getByToken(abs_failToken_,abs_failHandle_);
  const float &abs_fail = *abs_failHandle_;

  edm::Handle<std::vector<float>> chargeHandle_;
  iEvent.getByToken(chargeToken_,chargeHandle_);
  const std::vector<float> &charge = *chargeHandle_;

  std::cout << "Charge: " << std::endl;
  pcharge=0;
  for(auto c: charge){
    if (c > 0) pcharge=1;
  }
  std::cout << pcharge << std::endl;
  std::cout << "Endcharge" <<std::endl;

  fail = abs_fail;
  std::cout << "Failures: " << abs_fail << std::endl;

  std::map<float, GlobalPoint> gps_prop, gps_kf;
  std::map<float, float> map_xx_kf, map_xy_kf, map_yy_kf, map_xx_prop, map_xy_prop, map_yy_prop;

  std::ofstream kffile(outdir_ + "/Positions/kf_en"+energy_+"_eta_"+eta_+".csv", std::ios::app);

  std::cout << "# of KF points: " << kf.size() << std::endl;
  for(int i = 0;i<int(kf.size());i++){
    auto gp = kf[i];
    auto xx = xxkf[i];
    auto xy = xykf[i];
    auto yy = yykf[i];

    kffile << eventidx << "," << gp.x() << ","  << gp.y() << "," << gp.z() << "\n";
    gps_kf[gp.z()]=gp;
    map_xx_kf[gp.z()]=xx;
    map_xy_kf[gp.z()]=xy;
    map_yy_kf[gp.z()]=yy;
  }
  kffile.close();

  std::ofstream propfile(outdir_ + "/Positions/prop_en"+energy_+"_eta_"+eta_+".csv", std::ios::app);
  for(int i = 0;i<int(prop.size());i++){
    auto gp = prop[i];
    auto xx = xxprop[i];
    auto xy = xyprop[i];
    auto yy = yyprop[i];
    propfile << eventidx << ","  << gp.x() << ","  << gp.y() << "," << gp.z() << "\n";
    gps_prop[gp.z()] = gp;

    map_xx_prop[gp.z()]=xx;
    map_xy_prop[gp.z()]=xy;
    map_yy_prop[gp.z()]=yy;
  }
  
  propfile.close();

  // init vars
  const CaloGeometry &geom = iSetup.getData(caloGeomToken_);
  recHitTools_.setGeometry(geom);


  TString t_name;

  // Loop over Caloparticles 

  float edist = 0;
  float phidist = 0;
  float xdist = 0;
  float ydist = 0;
  float rdist = 0;
  float xpull = 0;
  float ypull = 0;

  std::ofstream simfile(outdir_+"/Positions/sim_en"+energy_+"_eta_"+eta_+".csv", std::ios::app);
  std::ofstream recfile(outdir_+"/Positions/rec_en"+energy_+"_eta_"+eta_+".csv", std::ios::app);
  std::ofstream efffile(outdir_+"/Efficiency/eff_en"+energy_+"_eta_"+eta_+".csv", std::ios::app);

  std::vector<DetId> tmprechits_; tmprechits_.clear();
  int rhits = 0;
  int smhits = 0;
  int eff_kf = 0;
  int eff_prop = 0;

  for (const auto& it_cp : cps) {
    // do something with track parameters, e.g, plot the charge.
    const CaloParticle& cp = ((it_cp)); 
    const SimClusterRefVector& simclusters = cp.simClusters();

    for (const auto& it_simc : simclusters){
      const SimCluster& simc = (*(it_simc));
      const auto& sc_hae = simc.hits_and_energies();

      for (const auto& it_sc_hae : sc_hae){

        // Characterize Hit

        DetId detid_ = (it_sc_hae.first);
        float hit_energy = it_sc_hae.second;
        std::map<DetId,const HGCRecHit *>::const_iterator itcheck = hitMap.find(detid_);
        unsigned int layer_ = recHitTools_.getLayerWithOffset(detid_); 

        std::string detector;
        std::string thickness;
        std::string tmp;

        if(recHitTools_.isSilicon(detid_)){
          detector = "Si";
          thickness = std::to_string(int(recHitTools_.getSiThickness(detid_))); 
          tmp = detector+" "+thickness;
          auto it = std::find(detectors.begin(), detectors.end(), tmp);
          if(it != detectors.end()){
            std::cout << "I'm alive" << std::endl;
            std::cout << tmp << std::endl;
          }
        }
        else{
          detector = "Sc";
          thickness = "None";
          tmp = "Sc";
          }

        // Position

        simfile << eventidx << ","  << recHitTools_.getPosition(detid_).x() << ","  << recHitTools_.getPosition(detid_).y() << "," << recHitTools_.getPosition(detid_).z() << "\n";

        // Hits & Energies


        layer_hits["Simhits"][""].front()->Fill(layer_,1);
        layer_hits["Simhits"][tmp].front()->Fill(layer_,1);
        hits["Simhits"][""].front()->Fill(1);
        hits["Simhits"][tmp].front()->Fill(1);
        layer_energies["Simhits"][""].front()->Fill(layer_, hit_energy);
        layer_energies["Simhits"][tmp].front()->Fill(layer_, hit_energy);
        energies["Simhits"][""].front()->Fill(hit_energy);
        energies["Simhits"][tmp].front()->Fill(hit_energy);
        single_layer_energies["Simhits"][tmp][layer_]->Fill(hit_energy);
        single_layer_energies["Simhits"][""][layer_]->Fill(hit_energy);
        single_layer_hits["Simhits"][tmp][layer_]->Fill(1);
        smhits=smhits+1;
        if (itcheck != hitMap.end()){
          rhits=rhits+1;
          layer_hits["Rechits"][""].front()->Fill(layer_,1);  
          layer_hits["Rechits"][tmp].front()->Fill(layer_,1);
          hits["Rechits"][""].front()->Fill(1);
          hits["Rechits"][tmp].front()->Fill(1);
          layer_energies["Rechits"][""].front()->Fill(layer_, hit_energy);
          layer_energies["Rechits"][tmp].front()->Fill(layer_, hit_energy);
          energies["Rechits"][""].front()->Fill(hit_energy);
          energies["Rechits"][tmp].front()->Fill(hit_energy);
          single_layer_energies["Rechits"][tmp][layer_]->Fill(hit_energy);
          single_layer_energies["Rechits"][""][layer_]->Fill(hit_energy);
          single_layer_hits["Rechits"][tmp][layer_]->Fill(1);
          //energies["Rechits"][detector].front()->Fill(hit_energy);
        }

        // KF & Propagator 
        for (const auto& pos: positions){
          auto &gps = (pos=="KF")? gps_kf:gps_prop;
          auto &map_xx = (pos=="KF")? map_xx_kf:map_xx_prop;
          auto &map_xy = (pos=="KF")? map_xx_kf:map_xx_prop;
          auto &map_yy = (pos=="KF")? map_xx_kf:map_xx_prop;
          auto &eff = (pos=="KF")? eff_kf:eff_prop;

          // Position

          auto gp = gps[recHitTools_.getPosition(detid_).z()];

          // Covariance Matrix

          float xx_disk = map_xx[recHitTools_.getPosition(detid_).z()];
          float xy_disk = map_xy[recHitTools_.getPosition(detid_).z()];
          float yy_disk = map_yy[recHitTools_.getPosition(detid_).z()];

          // Efficiency
          float e = 0;
          DetId closest_detid;
          if (detector == "Sc") closest_detid = static_cast<const HGCalGeometry*>(recHitTools_.getSubdetectorGeometry(detid_))->getClosestCell(gp);
          else closest_detid = static_cast<const HGCalGeometry*>(recHitTools_.getSubdetectorGeometry(detid_))->getClosestCellHex(gp, true);
          if (detid_ == closest_detid){
            std::map<DetId,const HGCRecHit *>::const_iterator echeck = hitMap.find(closest_detid);
            if (echeck != hitMap.end()){
              e = hitMap[closest_detid]->energy();
            }
            eff = eff+1;
            layer_hits[pos][""].front()->Fill(layer_,1);
            layer_hits[pos][tmp].front()->Fill(layer_,1);
            hits[pos][""].front()->Fill(1);
            hits[pos][tmp].front()->Fill(1);
            layer_energies[pos][""].front()->Fill(layer_, e);
            layer_energies[pos][tmp].front()->Fill(layer_, e);
            energies[pos][""].front()->Fill(e);
            energies[pos][tmp].front()->Fill(e);
            single_layer_energies[pos][""][layer_]->Fill(e);
            single_layer_energies[pos][tmp][layer_]->Fill(e);
            single_layer_hits[pos][tmp][layer_]->Fill(1);
          }

          // Residuals
          auto gp_det = recHitTools_.getPosition(detid_);
          if(gps.find(recHitTools_.getPosition(detid_).z())!=gps.end()){

              // Diff Eta Phi         
              edist = gp.eta() - gp_det.eta();
              phidist = gp.phi() - gp_det.phi();

              // Diff x y, pull
              xdist = gp.x() - gp_det.x();
              ydist = gp.y() - gp_det.y();

              xpull = xdist/xx_disk;
              ypull = ydist/yy_disk;

              // r Diff

              rdist = sqrt(gp.x()*gp.x() + gp.y()*gp.y()) - sqrt(gp_det.x()*gp_det.x() + gp_det.y()*gp_det.y());
              

              for(const auto& det : detectors){
                if(det == "" || det == detector || det.find(thickness)!=std::string::npos){
                  x_diff[pos][objects[0]][det].front()->Fill(xdist);
                  y_diff[pos][objects[0]][det].front()->Fill(ydist);
                  eta_diff[pos][objects[0]][det].front()->Fill(edist);
                  phi_diff[pos][objects[0]][det].front()->Fill(phidist); 
                  x_y_diff[pos][objects[0]][det].front()->Fill(xdist,ydist);
                  eta_phi_diff[pos][objects[0]][det].front()->Fill(edist,phidist);
                  r_diff[pos][objects[0]][det].front()->Fill(rdist);
                  x_pull[pos][objects[0]][det].front()->Fill(xpull);
                  y_pull[pos][objects[0]][det].front()->Fill(ypull);

                  layer_x_diff[pos][objects[0]][det].front()->Fill(layer_, xdist);
                  layer_y_diff[pos][objects[0]][det].front()->Fill(layer_, ydist);
                  layer_eta_diff[pos][objects[0]][det].front()->Fill(layer_, edist);
                  layer_phi_diff[pos][objects[0]][det].front()->Fill(layer_,phidist); 
                  layer_r_diff[pos][objects[0]][det].front()->Fill(layer_,rdist);

                  layer_abs_x_diff[pos][objects[0]][det].front()->Fill(layer_, std::abs(xdist));
                  layer_abs_y_diff[pos][objects[0]][det].front()->Fill(layer_, std::abs(ydist));
                  layer_abs_eta_diff[pos][objects[0]][det].front()->Fill(layer_, std::abs(edist));
                  layer_abs_phi_diff[pos][objects[0]][det].front()->Fill(layer_,std::abs(phidist)); 

                  layer_profile_abs_x_diff[pos][objects[0]][det].front()->Fill(layer_, std::abs(xdist),1);
                  layer_profile_abs_y_diff[pos][objects[0]][det].front()->Fill(layer_, std::abs(ydist),1);
                  layer_profile_abs_eta_diff[pos][objects[0]][det].front()->Fill(layer_, std::abs(edist),1);
                  layer_profile_abs_phi_diff[pos][objects[0]][det].front()->Fill(layer_,std::abs(phidist),1); 

                  layer_xydiff[pos][objects[0]][det][layer_]->Fill(xdist,ydist);

                  layer_eff[pos][objects[0]][det].front()->Fill(layer_, detid_ == closest_detid);
                }
              }
            
            
              if (itcheck != hitMap.end()){
                tmprechits_.push_back(detid_);
                recfile << eventidx << "," << recHitTools_.getPosition(detid_).x() << ","  << recHitTools_.getPosition(detid_).y() << "," << recHitTools_.getPosition(detid_).z() << "\n";

                for(const auto& det : detectors){
                  if(det == "" || det == detector || det.find(thickness)!=std::string::npos){

                    x_diff[pos][objects[1]][det].front()->Fill(xdist);
                    y_diff[pos][objects[1]][det].front()->Fill(ydist);
                    eta_diff[pos][objects[1]][det].front()->Fill(edist);
                    phi_diff[pos][objects[1]][det].front()->Fill(phidist); 
                    x_y_diff[pos][objects[1]][det].front()->Fill(xdist,ydist);
                    eta_phi_diff[pos][objects[1]][det].front()->Fill(edist,phidist);
                    r_diff[pos][objects[1]][det].front()->Fill(rdist);
                    x_pull[pos][objects[1]][det].front()->Fill(xpull);
                    y_pull[pos][objects[1]][det].front()->Fill(ypull);

                    layer_x_diff[pos][objects[1]][det].front()->Fill(layer_, xdist);
                    layer_y_diff[pos][objects[1]][det].front()->Fill(layer_, ydist);
                    layer_eta_diff[pos][objects[1]][det].front()->Fill(layer_, edist);
                    layer_phi_diff[pos][objects[1]][det].front()->Fill(layer_,phidist);
                    layer_r_diff[pos][objects[1]][det].front()->Fill(layer_,rdist);

                    layer_abs_x_diff[pos][objects[1]][det].front()->Fill(layer_, std::abs(xdist));
                    layer_abs_y_diff[pos][objects[1]][det].front()->Fill(layer_, std::abs(ydist));
                    layer_abs_eta_diff[pos][objects[1]][det].front()->Fill(layer_, std::abs(edist));
                    layer_abs_phi_diff[pos][objects[1]][det].front()->Fill(layer_,std::abs(phidist)); 
                
                    layer_profile_abs_x_diff[pos][objects[1]][det].front()->Fill(layer_, std::abs(xdist),1);
                    layer_profile_abs_y_diff[pos][objects[1]][det].front()->Fill(layer_, std::abs(ydist),1);
                    layer_profile_abs_eta_diff[pos][objects[1]][det].front()->Fill(layer_, std::abs(edist),1);
                    layer_profile_abs_phi_diff[pos][objects[1]][det].front()->Fill(layer_,std::abs(phidist),1); 

                    layer_xydiff[pos][objects[1]][det][layer_]->Fill(xdist,ydist);   

                    layer_eff[pos][objects[1]][det].front()->Fill(layer_, detid_ == closest_detid);
                  }
                }
              }
            }
          
          else{
            std::cout << "Incorrect z position " << recHitTools_.getPosition(detid_).z() << " at layer " << layer_<<std::endl;
          }
        }
      }
    }
  }
  std::cout << "Efficiency:" << eff_kf << std::endl;
  std::cout << "Hits:" << cps.size() << std::endl;
  efffile << eventidx << "," << eff_kf << "," << eff_prop << ","  << smhits << ","  << rhits << "\n";
  recfile.close();
  simfile.close();
  efffile.close();
  tree->Fill();

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void UpdatorStudies::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void UpdatorStudies::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void UpdatorStudies::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(UpdatorStudies);