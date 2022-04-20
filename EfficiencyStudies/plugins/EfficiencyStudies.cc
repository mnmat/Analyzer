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
#include <numeric>
#include <sstream>
#include <iomanip>

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


//ROOT includes
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"
#include <TFile.h>
#include <TROOT.h>
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
  void createStackPlot(TH1F*, TH1F*, TH1F*, std::string, std::string, std::string, std::string, std::string, std::string);
  void createTHSPlot(TCanvas*, std::vector<int>, std::vector<TH1F*>, TString, std::vector<TString>, std::vector<TString>, std::vector<float>, TString, TString);
  void createTH2Plot(TCanvas*, TH2F*, TString, std::vector<TString>, TString);
  void createTH1Plot(TCanvas*, TH1F*, TString, std::vector<TString>, TString);
  float getDr(float eta1, float phi1, float eta2, float phi2);
  std::vector<int> matchRecHit2CPRecHits(DetId detid_, std::vector<DetId> rechitdetid_);

  hgcal::RecHitTools recHitTools_;

  // ----------member data ---------------------------
  edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
  edm::EDGetTokenT<std::vector<CaloParticle>> caloParticlesToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsEEToken_; 
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsFHToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsBHToken_;
  edm::EDGetTokenT<reco::CaloClusterCollection> hgcalLayerClustersToken_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;

  TTree *tree = new TTree("tree","tree");
  // double trackPtMin_; 
  int skip_;
  std::string eta_;

  // Validate Particles

  TH1F *val_particle__sim_histo;
  TH1F *val_particle__rec_histo;
  TH1F *val_particle__lc_histo;

  // PID

  TH1F *pid_cp_histo;
  TH1F *pid_sim_histo;
  TH1F *pid_rec_histo;
  TH1F *pid_lc_histo;

  // Hit Plots

  int hit_sim;
  int hit_rec;
  int hit_lc;

  TH1F *hit_sim_histo;
  TH1F *hit_rec_histo;
  TH1F *hit_lc_histo;

  TH1F* hit_layer_sim_histo[47];
  TH1F* hit_layer_rec_histo[47];
  TH1F* hit_layer_lc_histo[47];

  TH1F* hit_layer_sim_si_histo[47];
  TH1F* hit_layer_rec_si_histo[47];
  TH1F* hit_layer_lc_si_histo[47];

  TH1F* hit_layer_sim_sc_histo[47];
  TH1F* hit_layer_rec_sc_histo[47];
  TH1F* hit_layer_lc_sc_histo[47];


  // Detector Hits Plots (Hits plotted over all layers/entire detector) 


  TH1F *det_sim_histo;
  TH1F *det_rec_histo;
  TH1F *det_lc_histo;

  TH1F *det_sim_si_histo;
  TH1F *det_rec_si_histo;
  TH1F *det_lc_si_histo;

  TH1F *det_sim_sc_histo;
  TH1F *det_rec_sc_histo;
  TH1F *det_lc_sc_histo;

  TH1F *det_bool_lc_histo;
  TH1F *det_bool_lc_si_histo;
  TH1F *det_bool_lc_sc_histo;

  TH1F *det_bool_rec_histo;
  TH1F *det_bool_rec_si_histo;
  TH1F *det_bool_rec_sc_histo;

  TH1F *det_bool_sim_histo;
  TH1F *det_bool_sim_si_histo;
  TH1F *det_bool_sim_sc_histo;

  // Connectivity plots

  TH1F* connectivity_lc_histo;
  TH1F* connectivity_rec_histo;
  TH1F* connectivity_sim_histo;

  TH1F* connectivity_lc_si_histo;
  TH1F* connectivity_rec_si_histo;
  TH1F* connectivity_sim_si_histo;

  TH1F* connectivity_lc_sc_histo;
  TH1F* connectivity_rec_sc_histo;
  TH1F* connectivity_sim_sc_histo;

  TH1F* connectivity_w_lc_histo;
  TH1F* connectivity_w_sim_histo;
  TH1F* connectivity_w_rec_histo;

  TH1F* connectivity_w_lc_si_histo;
  TH1F* connectivity_w_sim_si_histo;
  TH1F* connectivity_w_rec_si_histo;

  TH1F* connectivity_w_lc_sc_histo;
  TH1F* connectivity_w_sim_sc_histo;
  TH1F* connectivity_w_rec_sc_histo;

  TH1F* connectivity_max_lc_histo;
  TH1F* connectivity_max_sim_histo;
  TH1F* connectivity_max_rec_histo;

  TH1F* connectivity_max_lc_si_histo;
  TH1F* connectivity_max_sim_si_histo;
  TH1F* connectivity_max_rec_si_histo;

  TH1F* connectivity_max_lc_sc_histo;
  TH1F* connectivity_max_sim_sc_histo;
  TH1F* connectivity_max_rec_sc_histo;

  TH1F* connectivity_miss_lc_histo;
  TH1F* connectivity_miss_sim_histo;
  TH1F* connectivity_miss_rec_histo;

  TH1F* connectivity_miss_lc_si_histo;
  TH1F* connectivity_miss_sim_si_histo;
  TH1F* connectivity_miss_rec_si_histo;

  TH1F* connectivity_miss_lc_sc_histo;
  TH1F* connectivity_miss_sim_sc_histo;
  TH1F* connectivity_miss_rec_sc_histo;


  TH1F* connectivity_skip_lc_histo;
  TH1F* connectivity_skip_sim_histo;
  TH1F* connectivity_skip_rec_histo;

  TH1F* connectivity_skip_lc_si_histo;
  TH1F* connectivity_skip_sim_si_histo;
  TH1F* connectivity_skip_rec_si_histo;

  TH1F* connectivity_skip_lc_sc_histo;
  TH1F* connectivity_skip_sim_sc_histo;
  TH1F* connectivity_skip_rec_sc_histo;


  // Distance plots

  // TH1F* dist_eta_cp_simc;
  // TH1F* dist_phi_cp_simc;
  TH1F* dist_eta_cp_lc;
  TH1F* dist_phi_cp_lc;
  TH2F* dist_eta_phi_cp_lc;
  TH1F* dist_x_cp_lc;
  TH1F* dist_y_cp_lc;
  TH2F* dist_x_y_cp_lc;

  TH1F* dist_eta_cp_lc_si;
  TH1F* dist_phi_cp_lc_si;
  TH2F* dist_eta_phi_cp_lc_si;
  TH1F* dist_x_cp_lc_si;
  TH1F* dist_y_cp_lc_si;
  TH2F* dist_x_y_cp_lc_si;

  TH1F* dist_eta_cp_lc_sc;
  TH1F* dist_phi_cp_lc_sc;
  TH2F* dist_eta_phi_cp_lc_sc;
  TH1F* dist_x_cp_lc_sc;
  TH1F* dist_y_cp_lc_sc;
  TH2F* dist_x_y_cp_lc_sc;

  TH1F* eta_calo;
  TH1F* eta_sim;
  TH1F* eta_rec;
  TH1F* eta_lc;

  TH1F* dist_dr_cp_rec;
  TH1F* dist_eta_cp_rec;
  TH1F* dist_phi_cp_rec;
  TH2F* dist_eta_phi_cp_rec;
  TH1F* dist_x_cp_rec;
  TH1F* dist_y_cp_rec;
  TH2F* dist_x_y_cp_rec;

  TH1F* dist_dr_cp_rec_si;
  TH1F* dist_eta_cp_rec_si;
  TH1F* dist_phi_cp_rec_si;
  TH2F* dist_eta_phi_cp_rec_si;
  TH1F* dist_x_cp_rec_si;
  TH1F* dist_y_cp_rec_si;
  TH2F* dist_x_y_cp_rec_si;

  TH1F* dist_dr_cp_rec_sc;
  TH1F* dist_eta_cp_rec_sc;
  TH1F* dist_phi_cp_rec_sc;
  TH2F* dist_eta_phi_cp_rec_sc;
  TH1F* dist_x_cp_rec_sc;
  TH1F* dist_y_cp_rec_sc;
  TH2F* dist_x_y_cp_rec_sc;

  TH1F* dist_dr_cp_sim;
  TH1F* dist_eta_cp_sim;
  TH1F* dist_phi_cp_sim;
  TH2F* dist_eta_phi_cp_sim;
  TH1F* dist_x_cp_sim;
  TH1F* dist_y_cp_sim;
  TH2F* dist_x_y_cp_sim;

  TH1F* dist_dr_cp_sim_si;
  TH1F* dist_eta_cp_sim_si;
  TH1F* dist_phi_cp_sim_si;
  TH2F* dist_eta_phi_cp_sim_si;
  TH1F* dist_x_cp_sim_si;
  TH1F* dist_y_cp_sim_si;
  TH2F* dist_x_y_cp_sim_si;

  TH1F* dist_dr_cp_sim_sc;
  TH1F* dist_eta_cp_sim_sc;
  TH1F* dist_phi_cp_sim_sc;
  TH2F* dist_eta_phi_cp_sim_sc;
  TH1F* dist_x_cp_sim_sc;
  TH1F* dist_y_cp_sim_sc;
  TH2F* dist_x_y_cp_sim_sc;

  

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
      // trackPtMin_(iConfig.getParameter<double>("trackPtMin")), 
      caloGeomToken_(esConsumes<CaloGeometry, CaloGeometryRecord>()),
      skip_(iConfig.getParameter<int>("skip")),
      eta_(iConfig.getParameter<std::string>("eta")){



#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif

  //recHitTools_.reset(new hgcal::RecHitTools());
  //now do what ever initialization is needed


  usesResource("TFileService");
  edm::Service<TFileService> file;


  tree = file->make<TTree>("tree","hgc analyzer");

  // Branch (Review if this is necessary)
  
  /*
  
  tree->Branch("hit_sim", &hit_sim);
  tree->Branch("hit_rec", &hit_rec);
  tree->Branch("hit_lc", &hit_lc);

  for(int i = 0; i<47;i++){

    hit_layer_sim[i]=0;
    hit_layer_rec[i]=0;
    hit_layer_lc[i]=0;

    std::string sim_branch = "hit_layer_sim";
    std::string rec_branch = "hit_layer_rec";
    std::string lc_branch = "hit_layer_lc";

    std::stringstream nlayer;
    TString t_name;

    nlayer << i;
    sim_branch += nlayer.str();
    rec_branch += nlayer.str();
    lc_branch += nlayer.str();
    nlayer.str("");

    t_name = sim_branch;
    tree->Branch(t_name, &hit_layer_sim[i]);

    t_name = rec_branch;
    tree->Branch(t_name, &hit_layer_rec[i]);

    t_name = lc_branch;
    tree->Branch(t_name, &hit_layer_lc[i]);
  }

  */

  // Validation Plots (check Caloparticles, etc.)

  val_particle__sim_histo = new TH1F("Simclusters","Simclusters",5,0,5);
  val_particle__sim_histo->SetLineColor(kBlue);
//  val_particle__sim_histo->SetLineStyle(1);

  val_particle__rec_histo = new TH1F("Recclusters","Recclusters",500,0,500);
  val_particle__rec_histo->SetLineColor(kRed);
//  val_particle__rec_histo->SetLineStyle(2);

  val_particle__lc_histo = new TH1F("LCs","LCs",500,0,500);
  val_particle__lc_histo->SetLineColor(kGreen);
//  val_particle__lc_histo->SetLineStyle(7);


  // PID

  pid_cp_histo = new TH1F("PID CP","PID CP",40,-20,20);
  pid_sim_histo = new TH1F("PID Sim","PID Sim",40,-20,20);
  pid_rec_histo = new TH1F("PID Rec","PID Rec",40,-20,20);
  pid_lc_histo = new TH1F("PID LC","PID LC",40,-20,20);


  // Hit Plots

  hit_sim_histo = new TH1F("Simhits","Simhits",500,0,500);
  hit_sim_histo->SetLineColor(kBlue);
//  hit_sim_histo->SetLineStyle(1);

  hit_rec_histo = new TH1F("Rechits","Rechits",500,0,500);
  hit_rec_histo->SetLineColor(kRed);
//  hit_rec_histo->SetLineStyle(2);

  hit_lc_histo = new TH1F("LChits","LChits",500,0,500);
  hit_lc_histo->SetLineColor(kGreen);
//  hit_lc_histo->SetLineStyle(7);


  for(int i=0;i<47;i++){

    std::string sim_branch = "hit_layer_sim";
    std::string rec_branch = "hit_layer_rec";
    std::string lc_branch = "hit_layer_lc";

    std::stringstream nlayer;
    TString t_name;

    nlayer << i;
    sim_branch += nlayer.str();
    rec_branch += nlayer.str();
    lc_branch += nlayer.str();
    nlayer.str("");

    // Total hits

    t_name = sim_branch;
    hit_layer_sim_histo[i] = new TH1F(t_name,t_name,20,0,20);
    hit_layer_sim_histo[i]->SetLineColor(kBlue);
//    hit_layer_sim_histo[i]->SetLineStyle(1);


    t_name = rec_branch;
    hit_layer_rec_histo[i] = new TH1F(t_name,t_name,20,0,20);
    hit_layer_rec_histo[i]->SetLineColor(kRed);
//    hit_layer_rec_histo[i]->SetLineStyle(2);


    t_name = lc_branch;
    hit_layer_lc_histo[i] = new TH1F(t_name,t_name,20,0,20);
    hit_layer_lc_histo[i]->SetLineColor(kGreen);
//    hit_layer_lc_histo[i]->SetLineStyle(7);


    // Si hits

    t_name = sim_branch+"_si";
    hit_layer_sim_si_histo[i] = new TH1F(t_name,t_name,20,0,20);
    hit_layer_sim_si_histo[i]->SetLineColor(kBlue);
//    hit_layer_sim_si_histo[i]->SetLineStyle(1);


    t_name = rec_branch+"_si";
    hit_layer_rec_si_histo[i] = new TH1F(t_name,t_name,20,0,20);
    hit_layer_rec_si_histo[i]->SetLineColor(kRed);
//    hit_layer_rec_si_histo[i]->SetLineStyle(2);


    t_name = lc_branch+"_si";
    hit_layer_lc_si_histo[i] = new TH1F(t_name,t_name,20,0,20);
    hit_layer_lc_si_histo[i]->SetLineColor(kGreen);
//    hit_layer_lc_si_histo[i]->SetLineStyle(7);


    // Scintillator hits

    t_name = sim_branch+"_sc";
    hit_layer_sim_sc_histo[i] = new TH1F(t_name,t_name,20,0,20);
    hit_layer_sim_sc_histo[i]->SetLineColor(kBlue);
//    hit_layer_sim_sc_histo[i]->SetLineStyle(1);


    t_name = rec_branch+"_sc";
    hit_layer_rec_sc_histo[i] = new TH1F(t_name,t_name,20,0,20);
    hit_layer_rec_sc_histo[i]->SetLineColor(kRed);
//    hit_layer_rec_sc_histo[i]->SetLineStyle(2);


    t_name = lc_branch+"_sc";
    hit_layer_lc_sc_histo[i] = new TH1F(t_name,t_name,20,0,20);
    hit_layer_lc_sc_histo[i]->SetLineColor(kGreen);
//    hit_layer_lc_sc_histo[i]->SetLineStyle(7);

  }

  // Detector Hits Plots (Hits plotted over all layers/entire detector) 

  det_sim_histo = new TH1F("Simhits per Layer", "Simhits per Layer", 47, 0, 47);
  det_sim_si_histo = new TH1F("Simhits per Layer (Si)", "Simhits per Layer (Si)", 47, 0, 47);
  det_sim_sc_histo = new TH1F("Simhits per Layer (Scintillator)", "Simhits per Layer (Scintillator)", 47, 0, 47);

  det_rec_histo = new TH1F("Rechits per Layer", "Rechits per Layer", 47, 0, 47);
  det_rec_si_histo = new TH1F("Rechits per Layer (Si)", "Rechits per Layer (Si)", 47, 0, 47);
  det_rec_sc_histo = new TH1F("Rechits per Layer (Scintillator)", "Rechits per Layer (Scintillator)", 47, 0, 47);

  det_lc_histo = new TH1F("LChits per Layer", "LChits per Layer", 47, 0, 47);
  det_lc_si_histo = new TH1F("LChits per Layer (Si)", "LChits per Layer (Si)", 47, 0, 47);
  det_lc_sc_histo = new TH1F("LChits per Layer (Scintillator)", "LChits per Layer (Scintillator)", 47, 0, 47);

  det_bool_sim_histo = new TH1F("Sim layer total", "Sim layer total", 47,0,47);
  det_bool_sim_histo->SetLineColor(kBlue);
//  det_bool_sim_histo->SetLineStyle(1);

  det_bool_rec_histo = new TH1F("Rec layer total", "Rec layer total", 47,0,47);
  det_bool_rec_histo->SetLineColor(kRed);
//  det_bool_rec_histo->SetLineStyle(2);

  det_bool_lc_histo = new TH1F("LC layer total", "LC layer total", 47,0,47);
  det_bool_lc_histo->SetLineColor(kGreen);
//  det_bool_lc_histo->SetLineStyle(7);


  det_bool_sim_si_histo = new TH1F("Sim layer total (Si)", "Sim layer total (Si)", 47,0,47);
  det_bool_sim_si_histo->SetLineColor(kBlue);
//  det_bool_sim_si_histo->SetLineStyle(1);

  det_bool_rec_si_histo = new TH1F("Rec layer total (Si)", "Rec layer total (Si)", 47,0,47);
  det_bool_rec_si_histo->SetLineColor(kRed);
//  det_bool_rec_si_histo->SetLineStyle(2);

  det_bool_lc_si_histo = new TH1F("LC layer total (Si)", "LC layer total (Si)", 47,0,47);
  det_bool_lc_si_histo->SetLineColor(kGreen);
//  det_bool_lc_si_histo->SetLineStyle(7);


  det_bool_sim_sc_histo = new TH1F("Sim layer total (Sc)", "Sim layer total (Sc)", 47,0,47);
  det_bool_sim_sc_histo->SetLineColor(kBlue);
//  det_bool_sim_sc_histo->SetLineStyle(1);

  det_bool_rec_sc_histo = new TH1F("Rec layer total (Sc)", "Rec layer total (Sc)", 47,0,47);
  det_bool_rec_sc_histo->SetLineColor(kRed);
//  det_bool_rec_sc_histo->SetLineStyle(2);

  det_bool_lc_sc_histo = new TH1F("LC layer total (Sc)", "LC layer total (Sc)", 47,0,47);
  det_bool_lc_sc_histo->SetLineColor(kGreen);
//  det_bool_lc_sc_histo->SetLineStyle(7);


  // Connectivity Plots

  connectivity_lc_histo = new TH1F("LC Connectivity", "LC Connectivity", 48,0,48);
  connectivity_rec_histo = new TH1F("Rechit Connectivity", "Rechit Connectivity", 48,0,48);
  connectivity_sim_histo = new TH1F("Simhit Connectivity", "Simhit Connectivity", 48,0,48);

  connectivity_lc_si_histo = new TH1F("LC Connectivity (Si)", "LC Connectivity (Si)", 48,0,48);
  connectivity_rec_si_histo = new TH1F("Rechit Connectivity (Si)", "Rechit Connectivity (Si)", 48,0,48);
  connectivity_sim_si_histo = new TH1F("Simhit Connectivity (Si)", "Simhit Connectivity (Si)", 48,0,48);

  connectivity_lc_sc_histo = new TH1F("LC Connectivity (Sc)", "LC Connectivity (Sc)", 48,0,48);
  connectivity_rec_sc_histo = new TH1F("Rechit Connectivity (Sc)", "Rechit Connectivity (Sc)", 48,0,48);
  connectivity_sim_sc_histo = new TH1F("Simhit Connectivity (Sc)", "Simhit Connectivity (Sc)", 48,0,48);

  connectivity_w_lc_histo = new TH1F("LC Weighted Connectivity", "LC Weighted Connectivity", 48,0,48);
  connectivity_w_rec_histo = new TH1F("Rechit Weighted Connectivity", "Rechit Weighted Connectivity", 48,0,48);
  connectivity_w_sim_histo = new TH1F("Simhit Weighted Connectivity", "Simhit Weighted Connectivity", 48,0,48);

  connectivity_w_lc_si_histo = new TH1F("LC Weighted Connectivity (Si)", "LC Weighted Connectivity (Sc)", 48,0,48);
  connectivity_w_rec_si_histo = new TH1F("Rechit Weighted Connectivity (Si)", "Rechit Weighted Connectivity (Sc)", 48,0,48);
  connectivity_w_sim_si_histo = new TH1F("Simhit Weighted Connectivity (Si)", "Simhit Weighted Connectivity (Sc)", 48,0,48);

  connectivity_w_lc_sc_histo = new TH1F("LC Weighted Connectivity (Sc)", "LC Weighted Connectivity (Sc)", 48,0,48);
  connectivity_w_rec_sc_histo = new TH1F("Rechit Weighted Connectivity (Sc)", "Rechit Weighted Connectivity (Sc)", 48,0,48);
  connectivity_w_sim_sc_histo = new TH1F("Simhit Weighted Connectivity (Sc)", "Simhit Weighted Connectivity (Sc)", 48,0,48);

  connectivity_max_lc_histo = new TH1F("LC Max Connectivity", "LC Max Connectivity", 48,0,48);
  connectivity_max_rec_histo = new TH1F("Rechit Max Connectivity", "Rechit Max Connectivity", 48,0,48);
  connectivity_max_sim_histo = new TH1F("Simhit Max Connectivity", "Simhit Max Connectivity", 48,0,48);

  connectivity_max_lc_si_histo = new TH1F("LC Max Connectivity (Si)", "LC Max Connectivity (Si)", 48,0,48);
  connectivity_max_rec_si_histo = new TH1F("Rechit Max Connectivity (Si)", "Rechit Max Connectivity (Si)", 48,0,48);
  connectivity_max_sim_si_histo = new TH1F("Simhit Max Connectivity (Si)", "Simhit Max Connectivity (Si)", 48,0,48);

  connectivity_max_lc_sc_histo = new TH1F("LC Max Connectivity (Sc)", "LC Max Connectivity (Sc)", 48,0,48);
  connectivity_max_rec_sc_histo = new TH1F("Rechit Max Connectivity (Sc)", "Rechit Max Connectivity (Sc)", 48,0,48);
  connectivity_max_sim_sc_histo = new TH1F("Simhit Max Connectivity (Sc)", "Simhit Max Connectivity (Sc)", 48,0,48);

  connectivity_miss_lc_histo = new TH1F("LC Miss Connectivity", "LC Miss Connectivity", 48,0,48);
  connectivity_miss_rec_histo = new TH1F("Rechit Miss Connectivity", "Rechit Miss Connectivity", 48,0,48);
  connectivity_miss_sim_histo = new TH1F("Simhit Miss Connectivity", "Simhit Miss Connectivity", 48,0,48);

  connectivity_miss_lc_si_histo = new TH1F("LC Miss Connectivity (Si)", "LC Miss Connectivity (Si)", 48,0,48);
  connectivity_miss_rec_si_histo = new TH1F("Rechit Miss Connectivity (Si)", "Rechit Miss Connectivity (Si)", 48,0,48);
  connectivity_miss_sim_si_histo = new TH1F("Simhit Miss Connectivity (Si)", "Simhit Miss Connectivity (Si)", 48,0,48);

  connectivity_miss_lc_sc_histo = new TH1F("LC Miss Connectivity (Sc)", "LC Miss Connectivity (Sc)", 48,0,48);
  connectivity_miss_rec_sc_histo = new TH1F("Rechit Miss Connectivity (Sc)", "Rechit Miss Connectivity (Sc)", 48,0,48);
  connectivity_miss_sim_sc_histo = new TH1F("Simhit Miss Connectivity (Sc)", "Simhit Miss Connectivity (Sc)", 48,0,48);

  connectivity_skip_lc_histo = new TH1F("LC Skip Connectivity", "LC Skip Connectivity", 48,0,48);
  connectivity_skip_rec_histo = new TH1F("Rechit Skip Connectivity", "Rechit Skip Connectivity", 48,0,48);
  connectivity_skip_sim_histo = new TH1F("Simhit Skip Connectivity", "Simhit Skip Connectivity", 48,0,48);

  connectivity_skip_lc_si_histo = new TH1F("LC Skip Connectivity (Si)", "LC Skip Connectivity (Sc)", 48,0,48);
  connectivity_skip_rec_si_histo = new TH1F("Rechit Skip Connectivity (Si)", "Rechit Skip Connectivity (Si)", 48,0,48);
  connectivity_skip_sim_si_histo = new TH1F("Simhit Skip Connectivity (Si)", "Simhit Skip Connectivity (Si)", 48,0,48);

  connectivity_skip_lc_sc_histo = new TH1F("LC Skip Connectivity (Sc)", "LC Skip Connectivity (Sc)", 48,0,48);
  connectivity_skip_rec_sc_histo = new TH1F("Rechit Skip Connectivity (Sc)", "Rechit Skip Connectivity (Sc)", 48,0,48);
  connectivity_skip_sim_sc_histo = new TH1F("Simhit Skip Connectivity (Sc)", "Simhit Skip Connectivity (Sc)", 48,0,48);


  // Distance Plots


  // dist_eta_cp_simc = new TH1F("Eta Diff Simcluster - Caloparticles", "Eta Diff Simcluster - Caloparticles", 300,-0.15,0.15);
  // dist_phi_cp_simc = new TH1F("Phi Diff Simcluster - Caloparticles", "Phi Diff Simcluster - Caloparticles", 300,-0.15,0.15);


  eta_calo = new TH1F("Eta Caloparticles", "Eta Caloparticles", 300,0,3);
  eta_sim = new TH1F("Eta Simhits", "Eta Simhits", 300,0,3);
  eta_rec = new TH1F("Eta Rechits", "Eta Rechits", 300,0,3);
  eta_lc = new TH1F("Eta LCs", "Eta LCs", 300,0,3);

  dist_eta_cp_lc = new TH1F("Eta Diff Caloparticles - LC", "Eta Diff Caloparticles - LC", 300,-0.15,0.15);
  dist_phi_cp_lc = new TH1F("Phi Diff Caloparticles - LC", "Phi Diff Caloparticles - LC", 300,-0.15,0.15);
  dist_eta_phi_cp_lc = new TH2F("Eta Phi Diff Caloparticles - LC","Eta Phi Diff Caloparticles - LC",200,-0.2,0.2,200,-0.2,0.2);
  dist_x_cp_lc = new TH1F("x Diff Caloparticles - LC", "x Diff Caloparticles - LC", 200,-10,10);
  dist_y_cp_lc = new TH1F("y Diff Caloparticles - LC", "y Diff Caloparticles - LC", 200,-10,10);
  dist_x_y_cp_lc = new TH2F("x-y Diff Caloparticles - LC","x-y Diff Caloparticles - LC",200,-10,10,200,-10,10);

  dist_eta_cp_lc_si = new TH1F("Eta Diff Caloparticles - LC Si", "Eta Diff Caloparticles - LC Si", 300,-0.15,0.15);
  dist_phi_cp_lc_si = new TH1F("Phi Diff Caloparticles - LC Si", "Phi Diff Caloparticles - LC Si", 300,-0.15,0.15);
  dist_eta_phi_cp_lc_si = new TH2F("Eta Phi Diff Caloparticles - LC Si","Eta Phi Diff Caloparticles - LC Si",200,-0.2,0.2,200,-0.2,0.2);
  dist_x_cp_lc_si = new TH1F("x Diff Caloparticles - LC Si", "x Diff Caloparticles - LC Si", 200,-10,10);
  dist_y_cp_lc_si = new TH1F("y Diff Caloparticles - LC Si", "y Diff Caloparticles - LC Si", 200,-10,10);
  dist_x_y_cp_lc_si = new TH2F("x-y Diff Caloparticles - LC Si","x-y Diff Caloparticles - LC Si",200,-10,10,200,-10,10);

  dist_eta_cp_lc_sc = new TH1F("Eta Diff Caloparticles - LC Sc", "Eta Diff Caloparticles - LC Sc", 300,-0.15,0.15);
  dist_phi_cp_lc_sc = new TH1F("Phi Diff Caloparticles - LC Sc", "Phi Diff Caloparticles - LC Sc", 300,-0.15,0.15);
  dist_eta_phi_cp_lc_sc = new TH2F("Eta Phi Diff Caloparticles - LC Sc","Eta Phi Diff Caloparticles - LC Sc",200,-0.2,0.2,200,-0.2,0.2);
  dist_x_cp_lc_sc = new TH1F("x Diff Caloparticles - LC Sc", "x Diff Caloparticles - LC Sc", 200,-10,10);
  dist_y_cp_lc_sc = new TH1F("y Diff Caloparticles - LC Sc", "y Diff Caloparticles - LC Sc", 200,-10,10);
  dist_x_y_cp_lc_sc = new TH2F("x-y Diff Caloparticles - LC Sc","x-y Diff Caloparticles - LC Sc",200,-10,10,200,-10,10);

  dist_dr_cp_rec = new TH1F("Distance Caloparticles - Rechits","Distance Caloparticles - Rechits",100,0,1);
  dist_eta_cp_rec = new TH1F("Eta Diff Caloparticles - Rechits","Eta Diff Caloparticles - Rechits",200,-0.2,0.2);
  dist_phi_cp_rec = new TH1F("Phi Diff Caloparticles - Rechits","Phi Diff Caloparticles - Rechits",200,-0.2,0.2);
  dist_eta_phi_cp_rec = new TH2F("Eta Phi Diff Caloparticles - Rechits","Eta Phi Diff Caloparticles - Rechits",200,-0.2,0.2,200,-0.2,0.2);
  dist_x_cp_rec = new TH1F("x Diff Caloparticles - Rechits","x Diff Caloparticles - Rechits",200,-10,10);
  dist_y_cp_rec = new TH1F("y Diff Caloparticles - Rechits","y Diff Caloparticles - Rechits",200,-10,10);
  dist_x_y_cp_rec = new TH2F("x-y Diff Caloparticles - Rechits","x-y Diff Caloparticles - Rechits",200,-10,10,200,-10,10);

  dist_dr_cp_rec_si = new TH1F("Distance Caloparticles - Rechits Si","Distance Caloparticles - Rechits Si",100,0,1);
  dist_eta_cp_rec_si = new TH1F("Eta Diff Caloparticles - Rechits Si","Eta Diff Caloparticles - Rechits Si",200,-0.2,0.2);
  dist_phi_cp_rec_si = new TH1F("Phi Diff Caloparticles - Rechits Si","Phi Diff Caloparticles - Rechits Si",200,-0.2,0.2);
  dist_eta_phi_cp_rec_si = new TH2F("Eta Phi Diff Caloparticles - Rechits Si","Eta Phi Diff Caloparticles - Rechits Si",200,-0.2,0.2,200,-0.2,0.2);
  dist_x_cp_rec_si = new TH1F("x Diff Caloparticles - Rechits Si","x Diff Caloparticles - Rechits Si",200,-10,10);
  dist_y_cp_rec_si = new TH1F("y Diff Caloparticles - Rechits Si","y Diff Caloparticles - Rechits Si",200,-10,10);
  dist_x_y_cp_rec_si = new TH2F("x-y Diff Caloparticles - Rechits Si","x-y Diff Caloparticles - Rechits Si",200,-10,10,200,-10,10);

  dist_dr_cp_rec_sc = new TH1F("Distance Caloparticles - Rechits Sc","Distance Caloparticles - Rechits Sc",100,0,1);
  dist_eta_cp_rec_sc = new TH1F("Eta Diff Caloparticles - Rechits Sc","Eta Diff Caloparticles - Rechits Sc",200,-0.2,0.2);
  dist_phi_cp_rec_sc = new TH1F("Phi Diff Caloparticles - Rechits Sc","Phi Diff Caloparticles - Rechits Sc",200,-0.2,0.2);
  dist_eta_phi_cp_rec_sc = new TH2F("Eta Phi Diff Caloparticles - Rechits Sc","Eta Phi Diff Caloparticles - Rechits Sc",200,-0.2,0.2,200,-0.2,0.2);
  dist_x_cp_rec_sc = new TH1F("x Diff Caloparticles - Rechits Sc","x Diff Caloparticles - Rechits Sc",200,-10,10);
  dist_y_cp_rec_sc = new TH1F("y Diff Caloparticles - Rechits Sc","y Diff Caloparticles - Rechits Sc",200,-10,10);
  dist_x_y_cp_rec_sc = new TH2F("x-y Diff Caloparticles - Rechits Sc","x-y Diff Caloparticles - Rechits Sc",200,-10,10,200,-10,10);

  dist_dr_cp_sim = new TH1F("Distance Caloparticles - Simhits","Distance Caloparticles - Rechits",100,0,1);
  dist_eta_cp_sim = new TH1F("Eta Diff Caloparticles - Simhits","Eta Diff Caloparticles - Rechits",200,-0.2,0.2);
  dist_phi_cp_sim = new TH1F("Phi Diff Caloparticles - Simhits","Phi Diff Caloparticles - Rechits",200,-0.2,0.2);
  dist_eta_phi_cp_sim = new TH2F("Eta Phi Diff Caloparticles - Simhits","Eta Phi Diff Caloparticles - Rechits",200,-0.2,0.2,200,-0.2,0.2);
  dist_x_cp_sim = new TH1F("x Diff Caloparticles - Simhits","x Diff Caloparticles - Rechits",200,-10,10);
  dist_y_cp_sim = new TH1F("y Diff Caloparticles - Simhits","y DiffCaloparticles - Rechits",200,-10,10);
  dist_x_y_cp_sim = new TH2F("x-y Diff Caloparticles - Simhits","x-y Diff Caloparticles - Rechits",200,-10,10,200,-10,10);

  dist_dr_cp_sim_si = new TH1F("Distance Caloparticles - Simhits Si","Distance Caloparticles - Rechits Si",100,0,1);
  dist_eta_cp_sim_si = new TH1F("Eta Diff Caloparticles - Simhits Si","Eta Diff Caloparticles - Rechits Si",200,-0.2,0.2);
  dist_phi_cp_sim_si = new TH1F("Phi Diff Caloparticles - Simhits Si","Phi Diff Caloparticles - Rechits Si",200,-0.2,0.2);
  dist_eta_phi_cp_sim_si = new TH2F("Eta Phi Diff Caloparticles - Simhits Si","Eta Phi Diff Caloparticles - Rechits Si",200,-0.2,0.2,200,-0.2,0.2);
  dist_x_cp_sim_si = new TH1F("x Diff Caloparticles - Simhits Si","x Diff Caloparticles - Rechits Si",200,-10,10);
  dist_y_cp_sim_si = new TH1F("y Diff Caloparticles - Simhits Si","y DiffCaloparticles - Rechits Si",200,-10,10);
  dist_x_y_cp_sim_si = new TH2F("x-y Diff Caloparticles - Simhits Si","x-y Diff Caloparticles - Rechits Si",200,-10,10,200,-10,10);

  dist_dr_cp_sim_sc = new TH1F("Distance Caloparticles - Simhits Sc","Distance Caloparticles - Rechits Sc",100,0,1);
  dist_eta_cp_sim_sc = new TH1F("Eta Diff Caloparticles - Simhits Sc","Eta Diff Caloparticles - Rechits Sc",200,-0.2,0.2);
  dist_phi_cp_sim_sc = new TH1F("Phi Diff Caloparticles - Simhits Sc","Phi Diff Caloparticles - Rechits Sc",200,-0.2,0.2);
  dist_eta_phi_cp_sim_sc = new TH2F("Eta Phi Diff Caloparticles - Simhits Sc","Eta Phi Diff Caloparticles - Rechits Sc",200,-0.2,0.2,200,-0.2,0.2);
  dist_x_cp_sim_sc = new TH1F("x Diff Caloparticles - Simhits Sc","x Diff Caloparticles - Rechits Sc",200,-10,10);
  dist_y_cp_sim_sc = new TH1F("y Diff Caloparticles - Simhits Sc","y DiffCaloparticles - Rechits Sc",200,-10,10);
  dist_x_y_cp_sim_sc = new TH2F("x-y Diff Caloparticles - Simhits Sc","x-y Diff Caloparticles - Rechits Sc",200,-10,10,200,-10,10);

}

EfficiencyStudies::~EfficiencyStudies() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty


  //TLegend *leg = new TLegend(0.68,0.72,0.98,0.92); // Is dynamic memory allocation necessary?
  TCanvas *c1 = new TCanvas("c1","c1"); // Is dynamic memory allocation necessary?

  std::vector<int> colors = {1, 4, 2};
  TString folder = "/eos/user/m/mmatthew/www/Analyzer/"+eta_+"/";
  std::vector<float> pos = {0.70,0.75,0.9,0.9};

  // Validation Plots

  std::vector<TString> axes = {"# Particles", "# Events"};
  createTH1Plot(c1, val_particle__sim_histo, "SimParticles.png", axes, folder);

  axes = {"# Clusters","# Occurence"};
  createTH1Plot(c1, val_particle__lc_histo, "LayerClusters.png", axes, folder);


  // PID

  axes = {"# Clusters","PID"};
  createTH1Plot(c1, pid_cp_histo, "PID_CP.png", axes, folder);
  createTH1Plot(c1, pid_sim_histo, "PID_Sim.png", axes, folder);
  createTH1Plot(c1, pid_rec_histo, "PID_Rec.png", axes, folder);
  createTH1Plot(c1, pid_lc_histo, "PID_LC.png", axes, folder);

  // Hit Plots

  axes = {"Hits", "Events"};

  std::vector<TString> legend =  {"SimHits", "RecHits", "LCs"};
  std::vector<TH1F*> hists = {};
  TString tname;
  TString ttitle;

  for(int i=0; i<47; i++){
    // Total hits
    std::stringstream ss;
    ss << "Layer" <<std::setw(2) << std::setfill('2') << i <<".png";
    tname = ss.str();

    std::stringstream title;
    title << "Total Hits in Layer " << std::setfill('2') << i;
    ttitle = title.str();
    hists = {hit_layer_sim_histo[i], hit_layer_rec_histo[i], hit_layer_lc_histo[i]};
    createTHSPlot(c1, colors, hists, tname, axes, {"Simhits", "Rechits", "LCs"}, pos, folder, ttitle);

    // Si hits
    ss.str("");
    ss << "Layer" <<std::setw(2) << std::setfill('2') << i <<"_si.png";
    tname = ss.str();

    title.str("");
    title << "Si Hits in Layer " << std::setfill('2') << i;
    ttitle = title.str();
    hists = {hit_layer_sim_si_histo[i], hit_layer_rec_si_histo[i], hit_layer_lc_si_histo[i]};
    createTHSPlot(c1, colors, hists, tname, axes, {"Simhits", "Rechits", "LCs"}, pos, folder, ttitle);

    // Sc hits
    ss.str("");
    ss << "Layer" <<std::setw(2) << std::setfill('2') << i <<"_sc.png";
    tname = ss.str();

    title.str("");
    title << "Scintillator Hits in Layer " << std::setfill('2') << i;
    ttitle = title.str();
    hists = {hit_layer_sim_sc_histo[i], hit_layer_rec_sc_histo[i], hit_layer_lc_sc_histo[i]};
    createTHSPlot(c1, colors, hists, tname, axes, {"Simhits", "Rechits", "LCs"}, pos, folder, ttitle);
  } 

  createTH1Plot(c1, hit_sim_histo, "Simhits.png", axes, folder);
  createTH1Plot(c1, hit_rec_histo, "Rechits.png", axes, folder);
  createTH1Plot(c1, hit_lc_histo, "LChits.png", axes, folder);

  hists = {hit_sim_histo, hit_rec_histo, hit_lc_histo};
  createTHSPlot(c1, colors, hists, "Stackhits.png", axes, legend, pos, folder,"");


  // Det Hits Plots


  axes = {"Layer","# Clusters"};

  createTH1Plot(c1, det_bool_sim_histo,"LayerwiseSimClusters.png", axes, folder);
  createTH1Plot(c1, det_bool_lc_histo,"LayerwiseLCs.png", axes, folder);
  createTH1Plot(c1, det_bool_rec_histo,"LayerwiseRecClusters.png", axes, folder);

  axes = {"Layer","# Hits"};


  createTH1Plot(c1, det_sim_histo,"LayerwiseSimhits.png", axes, folder);
  createTH1Plot(c1, det_sim_si_histo,"LayerwiseSimhits_Si.png", axes, folder);
  createTH1Plot(c1, det_sim_sc_histo,"LayerwiseSimhits_Sc.png", axes, folder);

  hists = {det_sim_histo, det_sim_si_histo, det_sim_sc_histo};
  createTHSPlot(c1, colors, hists, "SimhitsStack.png", axes, {"Total", "Si", "Scintillator"}, pos, folder, "Test");

  createTH1Plot(c1, det_rec_histo,"LayerwiseRechits.png", axes, folder);
  createTH1Plot(c1, det_rec_si_histo,"LayerwiseRechits_Si.png", axes, folder);
  createTH1Plot(c1, det_rec_sc_histo,"LayerwiseRechits_Sc.png", axes, folder);

  hists = {det_rec_histo, det_rec_si_histo, det_rec_sc_histo};
  createTHSPlot(c1, colors, hists, "RechitsStack.png", axes, {"Total", "Si", "Scintillator"}, pos, folder, "Test");

  createTH1Plot(c1, det_lc_histo,"LayerwiseLChits.png", axes, folder);
  createTH1Plot(c1, det_lc_si_histo,"LayerwiseLChits_Si.png", axes, folder);
  createTH1Plot(c1, det_lc_sc_histo,"LayerwiseLChits_Sc.png", axes, folder);

  hists = {det_lc_histo, det_lc_si_histo, det_lc_sc_histo};
  createTHSPlot(c1, colors, hists, "LChitsStack.png", axes, {"Total", "Si", "Scintillator"}, pos, folder, "Test");

  // Connectivity Plots

  axes = {"Chain Length", "Occurence"};
  c1->SetLogy();

  createTH1Plot(c1, connectivity_lc_histo,"Connectivity_LC.png", axes, folder);
  createTH1Plot(c1, connectivity_rec_histo,"Connectivity_RecCluster.png", axes, folder);
  createTH1Plot(c1, connectivity_sim_histo,"Connectivity_SimCluster.png", axes, folder);

  connectivity_lc_histo ->SetStats(0);
  connectivity_rec_histo ->SetStats(0);
  connectivity_sim_histo ->SetStats(0);

  hists = {connectivity_sim_histo, connectivity_rec_histo, connectivity_lc_histo};
  createTHSPlot(c1, colors, hists, "Connectivity_Stack.png", axes, legend, {0.60,0.73,0.8,0.88}, folder, "Total Connectivity"); //{0.60,0.70,0.9,0.9}

  createTH1Plot(c1, connectivity_lc_si_histo,"Connectivity_LC_Si.png", axes, folder);
  createTH1Plot(c1, connectivity_rec_si_histo,"Connectivity_RecCluster_Si.png", axes, folder);
  createTH1Plot(c1, connectivity_sim_si_histo,"Connectivity_SimCluster_Si.png", axes, folder);

  connectivity_lc_si_histo ->SetStats(0);
  connectivity_rec_si_histo ->SetStats(0);
  connectivity_sim_si_histo ->SetStats(0);

  hists = {connectivity_sim_si_histo, connectivity_rec_si_histo, connectivity_lc_si_histo};
  createTHSPlot(c1, colors, hists, "Connectivity_Stack_Si.png", axes, legend, pos, folder, "Si Connectivity");

  createTH1Plot(c1, connectivity_lc_sc_histo,"Connectivity_LC_Sc.png", axes, folder);
  createTH1Plot(c1, connectivity_rec_sc_histo,"Connectivity_RecCluster_Sc.png", axes, folder);
  createTH1Plot(c1, connectivity_sim_sc_histo,"Connectivity_SimCluster_Sc.png", axes, folder);

  connectivity_lc_sc_histo ->SetStats(0);
  connectivity_rec_sc_histo ->SetStats(0);
  connectivity_sim_sc_histo ->SetStats(0);

  hists = {connectivity_sim_sc_histo, connectivity_rec_sc_histo, connectivity_lc_sc_histo};
  createTHSPlot(c1, colors, hists, "Connectivity_Stack_Scintillator.png", axes, legend, pos, folder, "Scintillator Connectivity");
  

  // Skip


  std::stringstream ss;
  ss << "Connectivity_Skip_" <<std::setw(1) << std::setfill('1') << skip_ << "_";
  tname = ss.str();

  createTH1Plot(c1, connectivity_skip_lc_histo,tname+"LC.png", axes, folder);
  createTH1Plot(c1, connectivity_skip_rec_histo,tname+"Rechit.png", axes, folder);
  createTH1Plot(c1, connectivity_skip_sim_histo,tname+"Simhit.png", axes, folder);

  connectivity_skip_lc_histo ->SetStats(0);
  connectivity_skip_rec_histo ->SetStats(0);
  connectivity_skip_sim_histo ->SetStats(0);

  hists = {connectivity_skip_sim_histo, connectivity_skip_rec_histo, connectivity_skip_lc_histo};
  createTHSPlot(c1, colors, hists, tname + "Stack.png", axes, legend, {0.60,0.73,0.8,0.88}, folder, "Skip Total Connectivity");

  createTH1Plot(c1, connectivity_skip_lc_si_histo,tname+"LC_Si.png", axes, folder);
  createTH1Plot(c1, connectivity_skip_rec_si_histo,tname+"RecCluster_Si.png", axes, folder);
  createTH1Plot(c1, connectivity_skip_sim_si_histo,tname+"SimCluster_Si.png", axes, folder);

  connectivity_skip_lc_si_histo ->SetStats(0);
  connectivity_skip_rec_si_histo ->SetStats(0);
  connectivity_skip_sim_si_histo ->SetStats(0);

  hists = {connectivity_skip_sim_si_histo, connectivity_skip_rec_si_histo, connectivity_skip_lc_si_histo};
  createTHSPlot(c1, colors, hists, tname+"Stack_Si.png", axes, legend, pos, folder, "Skip Si Connectivity");

  createTH1Plot(c1, connectivity_skip_lc_sc_histo,tname+"LC_Sc.png", axes, folder);
  createTH1Plot(c1, connectivity_skip_rec_sc_histo,tname+"RecCluster_Sc.png", axes, folder);
  createTH1Plot(c1, connectivity_skip_sim_sc_histo,tname+"SimCluster_Sc.png", axes, folder);

  connectivity_skip_lc_sc_histo ->SetStats(0);
  connectivity_skip_rec_sc_histo ->SetStats(0);
  connectivity_skip_sim_sc_histo ->SetStats(0);

  hists = {connectivity_skip_sim_sc_histo, connectivity_skip_rec_sc_histo, connectivity_skip_lc_sc_histo};
  createTHSPlot(c1, colors, hists, tname+"Stack_Scintillator.png", axes, legend, pos, folder, "Skip Scintillator Connectivity");

  c1->SetLogy(0);

  // Max

  createTH1Plot(c1, connectivity_max_lc_histo,"Connectivity_Max_LC.png", axes, folder);
  createTH1Plot(c1, connectivity_max_rec_histo,"Connectivity_Max_RecCluster.png", axes, folder);
  createTH1Plot(c1, connectivity_max_sim_histo,"Connectivity_Max_SimCluster.png", axes, folder);

  connectivity_max_lc_histo ->SetStats(0);
  connectivity_max_rec_histo ->SetStats(0);
  connectivity_max_sim_histo ->SetStats(0);

  hists = {connectivity_max_sim_histo, connectivity_max_rec_histo, connectivity_max_lc_histo};
  createTHSPlot(c1, colors, hists, "Connectivity_Max_Stack.png", axes, legend, {0.40,0.73,0.6,0.88}, folder, "Max Total Connectivity");

  createTH1Plot(c1, connectivity_max_lc_si_histo,"Connectivity_Max_LC_Si.png", axes, folder);
  createTH1Plot(c1, connectivity_max_rec_si_histo,"Connectivity_Max_RecCluster_Si.png", axes, folder);
  createTH1Plot(c1, connectivity_max_sim_si_histo,"Connectivity_Max_SimCluster_Si.png", axes, folder);

  connectivity_max_lc_si_histo ->SetStats(0);
  connectivity_max_rec_si_histo ->SetStats(0);
  connectivity_max_sim_si_histo ->SetStats(0);

  hists = {connectivity_max_sim_si_histo, connectivity_max_rec_si_histo, connectivity_max_lc_si_histo};
  createTHSPlot(c1, colors, hists, "Connectivity_Max_Stack_Si.png", axes, legend, pos, folder, "Max Si Connectivity");

  createTH1Plot(c1, connectivity_max_lc_sc_histo,"Connectivity_Max_LC_Sc.png", axes, folder);
  createTH1Plot(c1, connectivity_max_rec_sc_histo,"Connectivity_Max_RecCluster_Sc.png", axes, folder);
  createTH1Plot(c1, connectivity_max_sim_sc_histo,"Connectivity_Max_SimCluster_Sc.png", axes, folder);

  connectivity_max_lc_sc_histo ->SetStats(0);
  connectivity_max_rec_sc_histo ->SetStats(0);
  connectivity_max_sim_sc_histo ->SetStats(0);

  hists = {connectivity_max_sim_sc_histo, connectivity_max_rec_sc_histo, connectivity_max_lc_sc_histo};
  createTHSPlot(c1, colors, hists, "Connectivity_Max_Stack_Scintillator.png", axes, legend, pos, folder, "Max Scintillator Connectivity");

  // Miss

  createTH1Plot(c1, connectivity_miss_lc_histo,"Connectivity_Miss_LC.png", axes, folder);
  createTH1Plot(c1, connectivity_miss_rec_histo,"Connectivity_Miss_RecCluster.png", axes, folder);
  createTH1Plot(c1, connectivity_miss_sim_histo,"Connectivity_Miss_SimCluster.png", axes, folder);

  connectivity_miss_lc_histo ->SetStats(0);
  connectivity_miss_rec_histo ->SetStats(0);
  connectivity_miss_sim_histo ->SetStats(0);

  hists = {connectivity_miss_sim_histo, connectivity_miss_rec_histo, connectivity_miss_lc_histo};
  createTHSPlot(c1, colors, hists, "Connectivity_Miss_Stack.png", axes, legend, pos, folder, "Miss Total Connectivity");

  createTH1Plot(c1, connectivity_miss_lc_si_histo,"Connectivity_Miss_LC_Si.png", axes, folder);
  createTH1Plot(c1, connectivity_miss_rec_si_histo,"Connectivity_Miss_RecCluster_Si.png", axes, folder);
  createTH1Plot(c1, connectivity_miss_sim_si_histo,"Connectivity_Miss_SimCluster_Si.png", axes, folder);

  connectivity_miss_lc_si_histo ->SetStats(0);
  connectivity_miss_rec_si_histo ->SetStats(0);
  connectivity_miss_sim_si_histo ->SetStats(0);

  hists = {connectivity_miss_sim_si_histo, connectivity_miss_rec_si_histo, connectivity_miss_lc_si_histo};
  createTHSPlot(c1, colors, hists, "Connectivity_Miss_Stack_Si.png", axes, legend, pos, folder, "Miss Si Connectivity");

  createTH1Plot(c1, connectivity_miss_lc_sc_histo,"Connectivity_Miss_LC_Sc.png", axes, folder);
  createTH1Plot(c1, connectivity_miss_rec_sc_histo,"Connectivity_Miss_RecCluster_Sc.png", axes, folder);
  createTH1Plot(c1, connectivity_miss_sim_sc_histo,"Connectivity_Miss_SimCluster_Sc.png", axes, folder);

  connectivity_miss_lc_sc_histo ->SetStats(0);
  connectivity_miss_rec_sc_histo ->SetStats(0);
  connectivity_miss_sim_sc_histo ->SetStats(0);

  hists = {connectivity_miss_sim_sc_histo, connectivity_miss_rec_sc_histo, connectivity_miss_lc_sc_histo};
  createTHSPlot(c1, colors, hists, "Connectivity_Miss_Stack_Scintillator.png", axes, legend, pos, folder, "Miss Scintillator Connectivity");

  // Weighted

  axes = {"Layer", "Weighted Connectivity Score"};

  createTH1Plot(c1, connectivity_w_lc_histo,"WeightedConnectivity_LC.png",axes, folder);
  createTH1Plot(c1, connectivity_w_rec_histo,"WeightedConnectivity_RecCluster.png",axes, folder);
  createTH1Plot(c1, connectivity_w_sim_histo,"WeightedConnectivity_SimCluster.png",axes, folder);

  // Efficiency Plots

  axes = {"Layer", "Efficiency"};

  det_sim_histo->SetStats(0);
  TH1F *eff_sim_rec_histo = (TH1F*)det_rec_histo ->Clone("Efficiency Sim vs Rec");
  eff_sim_rec_histo->SetLineColor(kBlue);
  eff_sim_rec_histo->SetStats(0);
  TH1F *eff_sim_lc_histo = (TH1F*)det_lc_histo ->Clone("Efficiency Sim vs LC");
  eff_sim_lc_histo->SetLineColor(kRed);
  eff_sim_lc_histo->SetStats(0);

  eff_sim_rec_histo->Divide(det_sim_histo);
  eff_sim_lc_histo->Divide(det_sim_histo);

  hists = {eff_sim_rec_histo, eff_sim_lc_histo};
  createTHSPlot(c1, {4,2}, hists, "Efficiency_Plots.png", axes, {"RecHits", "LCs"}, {0.15,0.15,0.35,0.25}, folder, "Total Efficiency");

  TH1F *eff_sim_rec_si_histo = (TH1F*)det_rec_si_histo ->Clone("Efficiency Sim vs Rec (Si)");
  eff_sim_rec_si_histo->SetLineColor(kBlue);
  TH1F *eff_sim_lc_si_histo = (TH1F*)det_lc_si_histo ->Clone("Efficiency Sim vs LC (Si)");
  eff_sim_lc_si_histo->SetLineColor(kRed);

  eff_sim_rec_si_histo->Divide(det_sim_si_histo);
  eff_sim_lc_si_histo->Divide(det_sim_si_histo);

  hists = {eff_sim_rec_si_histo, eff_sim_lc_si_histo};
  createTHSPlot(c1, {4,2}, hists, "Efficiency_Plots_Si.png", axes, {"RecHits", "Lcs"},{0.15,0.15,0.35,0.25}, folder, "Si Efficiency");

  TH1F *eff_sim_rec_sc_histo = (TH1F*)det_rec_sc_histo ->Clone("Efficiency Sim vs Rec (Scintillator)");
  eff_sim_rec_sc_histo->SetLineColor(kRed);
  TH1F *eff_sim_lc_sc_histo = (TH1F*)det_lc_sc_histo ->Clone("Efficiency Sim vs LC (Scintillator)");
  eff_sim_lc_sc_histo->SetLineColor(kBlue);

  eff_sim_rec_sc_histo->Divide(det_sim_sc_histo);
  eff_sim_lc_sc_histo->Divide(det_sim_sc_histo);

  hists = {eff_sim_rec_sc_histo, eff_sim_lc_sc_histo};
  createTHSPlot(c1, {4,2}, hists, "Efficiency_Plots_Scintillator.png", axes, {"RecHits", "LCs"},{0.15,0.15,0.35,0.25}, folder, "Scintillator Efficiency");

  // Distance Plots

  // createTH1Plot(c1, dist_eta_cp_simc, "EtaSimcCalo.png", axes, folder);
  // createTH1Plot(c1, dist_phi_cp_simc, "PhiSimcCalo.png", axes, folder);

  createTH1Plot(c1, eta_calo, "Eta_Calo.png", axes, folder);
  createTH1Plot(c1, eta_sim, "Eta_Sim.png", axes, folder);
  createTH1Plot(c1, eta_rec, "Eta_Rec.png", axes, folder);
  createTH1Plot(c1, eta_lc, "Eta_LC.png", axes, folder);


  axes = {"#Delta R","# Occurence"};
  createTH1Plot(c1, dist_dr_cp_rec, "Dist_Sim_CP.png", axes, folder);

  axes = {"#Delta #eta", "# Occurence"};
  createTH1Plot(c1, dist_eta_cp_sim,"Diff_Eta_CP_Sim.png", axes, folder);
  createTH1Plot(c1, dist_eta_cp_rec,"Diff_Eta_CP_Rec.png", axes, folder);
  createTH1Plot(c1, dist_eta_cp_lc, "Diff_Eta_CP_LC.png", axes, folder);

  createTH1Plot(c1, dist_eta_cp_sim_si,"Diff_Eta_CP_Sim_Si.png", axes, folder);
  createTH1Plot(c1, dist_eta_cp_rec_si,"Diff_Eta_CP_Rec_Si.png", axes, folder);
  createTH1Plot(c1, dist_eta_cp_lc_si, "Diff_Eta_CP_LC_Si.png", axes, folder);

  createTH1Plot(c1, dist_eta_cp_sim_sc,"Diff_Eta_CP_Sim_Sc.png", axes, folder);
  createTH1Plot(c1, dist_eta_cp_rec_sc,"Diff_Eta_CP_Rec_Sc.png", axes, folder);
  createTH1Plot(c1, dist_eta_cp_lc_sc, "Diff_Eta_CP_LC_Sc.png", axes, folder);

  axes = {"#Delta #phi", "# Occurence"};
  createTH1Plot(c1, dist_phi_cp_sim,"Diff_Phi_CP_Sim.png", axes, folder);
  createTH1Plot(c1, dist_phi_cp_rec,"Diff_Phi_CP_Rec.png", axes, folder);
  createTH1Plot(c1, dist_phi_cp_lc, "Diff_Phi_CP_LC.png", axes, folder);

  createTH1Plot(c1, dist_phi_cp_sim_si,"Diff_Phi_CP_Sim_Si.png", axes, folder);
  createTH1Plot(c1, dist_phi_cp_rec_si,"Diff_Phi_CP_Rec_Si.png", axes, folder);
  createTH1Plot(c1, dist_phi_cp_lc_si, "Diff_Phi_CP_LC_Si.png", axes, folder);

  createTH1Plot(c1, dist_phi_cp_sim_sc,"Diff_Phi_CP_Sim_Sc.png", axes, folder);
  createTH1Plot(c1, dist_phi_cp_rec_sc,"Diff_Phi_CP_Rec_Sc.png", axes, folder);
  createTH1Plot(c1, dist_phi_cp_lc_sc, "Diff_Phi_CP_LC_Sc.png", axes, folder);

  axes = {"#Delta #eta", "#Delta #phi"};
  createTH2Plot(c1,dist_eta_phi_cp_sim,"Diff_Eta_Phi_CP_Sim.png", axes, folder);
  createTH2Plot(c1,dist_eta_phi_cp_rec,"Diff_Eta_Phi_CP_Rec.png", axes, folder);
  createTH2Plot(c1,dist_eta_phi_cp_lc,"Diff_Eta_Phi_CP_LC.png", axes, folder);

  createTH2Plot(c1,dist_eta_phi_cp_sim_si,"Diff_Eta_Phi_CP_Sim_Si.png", axes, folder);
  createTH2Plot(c1,dist_eta_phi_cp_rec_si,"Diff_Eta_Phi_CP_Rec_Si.png", axes, folder);
  createTH2Plot(c1,dist_eta_phi_cp_lc_si,"Diff_Eta_Phi_CP_LC_Si.png", axes, folder);

  createTH2Plot(c1,dist_eta_phi_cp_sim,"Diff_Eta_Phi_CP_Sim_Sc.png", axes, folder);
  createTH2Plot(c1,dist_eta_phi_cp_rec,"Diff_Eta_Phi_CP_Rec_Sc.png", axes, folder);
  createTH2Plot(c1,dist_eta_phi_cp_lc,"Diff_Eta_Phi_CP_LC_Sc.png", axes, folder);

  axes = {"#Delta x", "# Occurence"};
  createTH1Plot(c1, dist_x_cp_sim,"Diff_x_CP_Sim.png", axes, folder);
  createTH1Plot(c1, dist_x_cp_rec,"Diff_x_CP_Rec.png", axes, folder);
  createTH1Plot(c1, dist_x_cp_lc, "Diff_x_CP_LC.png", axes, folder);

  createTH1Plot(c1, dist_x_cp_sim_si,"Diff_x_CP_Sim_Si.png", axes, folder);
  createTH1Plot(c1, dist_x_cp_rec_si,"Diff_x_CP_Rec_Si.png", axes, folder);
  createTH1Plot(c1, dist_x_cp_lc_si, "Diff_x_CP_LC_Si.png", axes, folder);

  createTH1Plot(c1, dist_x_cp_sim_sc,"Diff_x_CP_Sim_Sc.png", axes, folder);
  createTH1Plot(c1, dist_x_cp_rec_sc,"Diff_x_CP_Rec_Sc.png", axes, folder);
  createTH1Plot(c1, dist_x_cp_lc_sc, "Diff_x_CP_LC_Sc.png", axes, folder);

  axes = {"#Delta y", "# Occurence"};
  createTH1Plot(c1, dist_y_cp_sim,"Diff_y_CP_Sim.png", axes, folder);
  createTH1Plot(c1, dist_y_cp_rec,"Diff_y_CP_Rec.png", axes, folder);
  createTH1Plot(c1, dist_y_cp_lc, "Diff_y_CP_LC.png", axes, folder);

  createTH1Plot(c1, dist_y_cp_sim_si,"Diff_y_CP_Sim_Si.png", axes, folder);
  createTH1Plot(c1, dist_y_cp_rec_si,"Diff_y_CP_Rec_Si.png", axes, folder);
  createTH1Plot(c1, dist_y_cp_lc_si, "Diff_y_CP_LC_Si.png", axes, folder);

  createTH1Plot(c1, dist_y_cp_sim_sc,"Diff_y_CP_Sim_Sc.png", axes, folder);
  createTH1Plot(c1, dist_y_cp_rec_sc,"Diff_y_CP_Rec_Sc.png", axes, folder);
  createTH1Plot(c1, dist_y_cp_lc_sc, "Diff_y_CP_LC_Sc.png", axes, folder);

  axes = {"#Delta x", "#Delta y"};
  createTH2Plot(c1,dist_x_y_cp_sim,"Diff_xy_CP_Sim.png", axes, folder);
  createTH2Plot(c1,dist_x_y_cp_rec,"Diff_xy_CP_Rec.png", axes, folder);
  createTH2Plot(c1,dist_x_y_cp_lc,"Diff_xy_CP_LC.png", axes, folder);

  createTH2Plot(c1,dist_x_y_cp_sim_si,"Diff_xy_CP_Sim_Si.png", axes, folder);
  createTH2Plot(c1,dist_x_y_cp_rec_si,"Diff_xy_CP_Rec_Si.png", axes, folder);
  createTH2Plot(c1,dist_x_y_cp_lc_si,"Diff_xy_CP_LC_Si.png", axes, folder);

  createTH2Plot(c1,dist_x_y_cp_sim_sc,"Diff_xy_CP_Sim_Sc.png", axes, folder);
  createTH2Plot(c1,dist_x_y_cp_rec_sc,"Diff_xy_CP_Rec_Sc.png", axes, folder);
  createTH2Plot(c1,dist_x_y_cp_lc_sc,"Diff_xy_CP_LC_Sc.png", axes, folder);
}



//
// member functions
//

// ------------ method called for each event  ------------

void EfficiencyStudies::createTH1Plot(TCanvas *c, TH1F *hs, TString name, std::vector<TString> axis, TString folder){

  gSystem->cd(folder);
  c->cd();
  hs->GetXaxis()->SetTitle(axis[0]);
  hs->GetYaxis()->SetTitle(axis[1]);
  hs->Write();
  hs->Draw();
  c->SaveAs(name);
}

void EfficiencyStudies::createTH2Plot(TCanvas *c, TH2F *hs, TString name, std::vector<TString> axis, TString folder){

  gSystem->cd(folder);
  c->cd();
  hs->GetXaxis()->SetTitle(axis[0]);
  hs->GetYaxis()->SetTitle(axis[1]);
  hs->Write();
  hs->Draw();
  c->SaveAs(name);
}

void EfficiencyStudies::createTHSPlot(TCanvas *c, std::vector<int> colors, std::vector<TH1F*> hist, TString name, std::vector<TString> axis, std::vector<TString> lname, std::vector<float> lpos, TString folder,TString title){

  THStack hs = THStack(title,title);
  TLegend l = TLegend(lpos[0],lpos[1],lpos[2],lpos[3]);
  std::vector<int> styles = {1,2,7};
  std::vector<int> widths = {3,2,1};

  gSystem->cd(folder); // move outside?
  c->cd(); // move outside?

  for(int i=0; i<int(colors.size());i++){
    hist[i]->SetLineColor(colors[i]);
    hist[i]->SetLineStyle(styles[i]);
    hist[i]->SetLineWidth(widths[i]);
    hist[i]->Write();
    hs.Add(hist[i]);
    l.AddEntry(hist[i],lname[i],"l");
  }

  hs.Draw("nostack");
  gPad->Modified();
  gPad->Update();

  hs.GetXaxis()->SetTitle(axis[0]);
  hs.GetYaxis()->SetTitle(axis[1]);

  gPad->Modified();
  gPad->Update();
  
  l.Draw();
  c->SaveAs(name);
}


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

std::vector<int> EfficiencyStudies::matchRecHit2CPRecHits(DetId detid_, std::vector<DetId> rechitdetid_) {
  std::vector<int> matchedIdxs; matchedIdxs.clear();
  for (unsigned int i0=0; i0<rechitdetid_.size(); ++i0) {
    if (detid_ == rechitdetid_[i0]) { matchedIdxs.push_back(i0); }
  }
  return matchedIdxs;
} // end of matchRecHit2CPRecHits


float EfficiencyStudies::getDr(float eta1, float phi1, float eta2, float phi2) { 
  return sqrt( (eta1-eta2)*(eta1-eta2) + (phi1-phi2)*(phi1-phi2) ); 
}


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
  recHitTools_.setGeometry(geom);

  int nTrack = 0;
  for (const auto& track : iEvent.get(tracksToken_)) {
    // do something with track parameters, e.g, plot the charge.
    // int charge = track.charge();
    // if(track.pt() < trackPtMin_) continue;
    nTrack++;
  }

  int nCalo = 0;
  int nSim = 0;
  int nLC = 0;
  int nSimhits = 0;
  int nRechits = 0;
  int nLChits = 0;

  std::vector<int> hit_layer_sim(47,0);
  std::vector<int> hit_layer_sim_si(47,0);
  std::vector<int> hit_layer_sim_sc(47,0);
  std::vector<int> hit_layer_rec(47,0);
  std::vector<int> hit_layer_rec_si(47,0);
  std::vector<int> hit_layer_rec_sc(47,0);
  std::vector<int> hit_layer_lc(47,0);
  std::vector<int> hit_layer_lc_si(47,0);
  std::vector<int> hit_layer_lc_sc(47,0);

  /*

  for(int i=0;i<47;i++){
    hit_layer_sim[i]=0;
    hit_layer_sim_si[i]=0;
    hit_layer_sim_sc[i]=0;

    hit_layer_rec[i]=0;
    hit_layer_rec_si[i]=0;
    hit_layer_rec_sc[i]=0;

    hit_layer_lc[i]=0;
    hit_layer_lc_si[i]=0;
    hit_layer_lc_sc[i]=0;
    }

  */

  std::vector<double> temp_simhits(47,0);
  std::vector<double> temp_rechits(47,0);
  std::vector<double> temp_lchits(47,0);

  float cp_eta = 0;
  float cp_phi = 0;


  // Loop over Caloparticles 


  std::vector<DetId> tmprechits_; tmprechits_.clear();
  for (const auto& it_cp : cps) {
    // do something with track parameters, e.g, plot the charge.
    // int charge = track.charge();
    nCalo++;
    const CaloParticle& cp = ((it_cp)); 

    cp_eta = it_cp.eta();
    cp_phi = it_cp.phi();

    eta_calo -> Fill(cp_eta); // Investigation of Scintillator-Eta 2.9 discrepancy
    pid_cp_histo -> Fill(cp.pdgId());

    const SimClusterRefVector& simclusters = cp.simClusters();
    // std::cout<<"CP number of Simhits:" << simclusters.numberOfSimHits() << std::endl;    

    for (const auto& it_simc : simclusters){

      nSim++;

      const SimCluster& simc = (*(it_simc));

      const std::vector<SimTrack>& simt = simc.g4Tracks();
      for (const auto& it_simt : simt){
        std::cout << it_simt <<endl;
        pid_sim_histo -> Fill(it_simt.type());
      }


      // Investigating simcluster position
      
      // dist_eta_cp_simc->Fill(cp_eta - simc.eta());
      // dist_phi_cp_simc->Fill(cp_phi - simc.phi());

      const auto& sc_haf = simc.hits_and_fractions();
      //std::cout<<"CP number of Simhits:" << simc.numberOfSimHits() << std::endl;
      std::cout<<"CP number of Rechits:" << simc.numberOfRecHits() << std::endl;

      for (const auto& it_sc_haf : sc_haf){
        nSimhits++;
        DetId detid_ = (it_sc_haf.first);
        std::map<DetId,const HGCRecHit *>::const_iterator itcheck = hitMap.find(detid_);
        unsigned int layer_ = recHitTools_.getLayerWithOffset(detid_);   


        // Investigation of Scintillator-Eta 2.9 discrepancy


        float edist = cp_eta - recHitTools_.getPosition(detid_).eta();
        float phidist = cp_phi - recHitTools_.getPosition(detid_).phi();
        float drdist = getDr(recHitTools_.getPosition(detid_).eta(),recHitTools_.getPosition(detid_).phi(),cp_eta,cp_phi);

        float cpPositionAtSimHit_x = -1.*recHitTools_.getPosition(detid_).z()*TMath::Tan(2.*(TMath::ATan(TMath::Exp(cp_eta))))*TMath::Cos(cp_phi);
        float cpPositionAtSimHit_y = -1.*recHitTools_.getPosition(detid_).z()*TMath::Tan(2.*(TMath::ATan(TMath::Exp(cp_eta))))*TMath::Sin(cp_phi);

        float xdist = cpPositionAtSimHit_x - recHitTools_.getPosition(detid_).x();
        float ydist = cpPositionAtSimHit_y - recHitTools_.getPosition(detid_).y();

        eta_sim -> Fill(recHitTools_.getPosition(detid_).eta());
        
        dist_dr_cp_sim->Fill(drdist);
        dist_eta_cp_sim->Fill(edist);
        dist_phi_cp_sim->Fill(phidist);
        dist_eta_phi_cp_sim->Fill(edist,phidist); 

        dist_x_cp_sim->Fill(xdist);
        dist_y_cp_sim->Fill(ydist);
        dist_x_y_cp_sim->Fill(xdist,ydist);


        // End of investigation
      
        if(recHitTools_.isSilicon(detid_)){
          hit_layer_sim_si[layer_-1]++;

          dist_dr_cp_sim_si->Fill(drdist);
          dist_eta_cp_sim_si->Fill(edist);
          dist_phi_cp_sim_si->Fill(phidist);
          dist_eta_phi_cp_sim_si->Fill(edist,phidist); 

          dist_x_cp_sim_si->Fill(xdist);
          dist_y_cp_sim_si->Fill(ydist);
          dist_x_y_cp_sim_si->Fill(xdist,ydist);
          
          }
        if(recHitTools_.isScintillator(detid_)){
          hit_layer_sim_sc[layer_-1]++;
          std::cout<<"Scintillator hit"<<std::endl;

          dist_dr_cp_sim_sc->Fill(drdist);
          dist_eta_cp_sim_sc->Fill(edist);
          dist_phi_cp_sim_sc->Fill(phidist);
          dist_eta_phi_cp_sim_sc->Fill(edist,phidist); 

          dist_x_cp_sim_sc->Fill(xdist);
          dist_y_cp_sim_sc->Fill(ydist);
          dist_x_y_cp_sim_sc->Fill(xdist,ydist);
        }
    
	      hit_layer_sim[layer_-1]++;

        if (itcheck != hitMap.end()){
          tmprechits_.push_back(detid_);
          nRechits++;

          // Not necessary? Detid of Simhit and Rechit the same.

		      edist = cp_eta - recHitTools_.getPosition(detid_).eta();
		      phidist = cp_phi - recHitTools_.getPosition(detid_).phi();
		      drdist = getDr(recHitTools_.getPosition(detid_).eta(),recHitTools_.getPosition(detid_).phi(),cp_eta,cp_phi);

          GlobalPoint rhPosition_ = recHitTools_.getPosition(detid_);
          float rhPosition_x = rhPosition_.x();
          float rhPosition_y = rhPosition_.y();
        
          float cpPositionAtRecHit_x = -1.*rhPosition_.z()*TMath::Tan(2.*(TMath::ATan(TMath::Exp(cp_eta))))*TMath::Cos(cp_phi);
          float cpPositionAtRecHit_y = -1.*rhPosition_.z()*TMath::Tan(2.*(TMath::ATan(TMath::Exp(cp_eta))))*TMath::Sin(cp_phi);

          float xdist = cpPositionAtSimHit_x - recHitTools_.getPosition(detid_).x();
          float ydist = cpPositionAtSimHit_y - recHitTools_.getPosition(detid_).y();

          //float drRecHitCP = getDr(rhPosition_x,rhPosition_y,cpPositionAtRecHit_x,cpPositionAtRecHit_y);


          if(recHitTools_.isSilicon(detid_)){hit_layer_rec_si[layer_-1]++;
            //dist_dr_cp_rec_si->Fill(drdist);
            dist_eta_cp_rec_si->Fill(edist);
            dist_phi_cp_rec_si->Fill(phidist);
            dist_eta_phi_cp_rec_si->Fill(edist,phidist); 

            dist_x_cp_rec_si->Fill(xdist);
            dist_y_cp_rec_si->Fill(ydist);
            dist_x_y_cp_rec_si->Fill(xdist,ydist);
          }
          if(recHitTools_.isScintillator(detid_)){
            hit_layer_rec_sc[layer_-1]++;

            //dist_dr_cp_rec_sc->Fill(drdist);
            dist_eta_cp_rec_sc->Fill(edist);
            dist_phi_cp_rec_sc->Fill(phidist);
            dist_eta_phi_cp_rec_sc->Fill(edist,phidist); 

            dist_x_cp_rec_sc->Fill(xdist);
            dist_y_cp_rec_sc->Fill(ydist);
            dist_x_y_cp_rec_sc->Fill(xdist,ydist);
          }


          hit_layer_rec[layer_-1]++;


		      dist_dr_cp_rec->Fill(drdist);
		      dist_eta_cp_rec->Fill(edist);
		      dist_phi_cp_rec->Fill(phidist);
		      dist_eta_phi_cp_rec->Fill(edist,phidist); 


          dist_x_cp_rec -> Fill(cpPositionAtRecHit_x-rhPosition_x);
          dist_y_cp_rec -> Fill(cpPositionAtRecHit_y-rhPosition_y);
          dist_x_y_cp_rec -> Fill(cpPositionAtRecHit_x-rhPosition_x, cpPositionAtRecHit_y-rhPosition_y);

          //dist_dr_cp_rec->Fill(drRecHitCP);
        }
      }
    }
  }
  

  // Loop over LCs  

  
  for (const auto& it_lc : lcs) {
    nLC++;
    const std::vector<std::pair<DetId, float>> &hf = it_lc.hitsAndFractions();


    // loop over the rechits of this specific layer cluster
    for (unsigned int j = 0; j < hf.size(); j++) { 
      DetId detid_ = hf[j].first;
      std::vector<int> idx_= matchRecHit2CPRecHits(detid_, tmprechits_);

      for (unsigned int i0=0; i0<idx_.size(); i0++) {
        nLChits++;
        unsigned int layer_ = recHitTools_.getLayerWithOffset(tmprechits_[idx_[i0]]); 

	      float edist = cp_eta - recHitTools_.getPosition(detid_).eta();
	      float phidist = cp_phi - recHitTools_.getPosition(detid_).phi();
	      float drdist = getDr(recHitTools_.getPosition(detid_).eta(),recHitTools_.getPosition(detid_).phi(),cp_eta,cp_phi);

        float cpPositionAtLC_x = -1.*it_lc.z()*TMath::Tan(2.*(TMath::ATan(TMath::Exp(cp_eta))))*TMath::Cos(cp_phi);
        float cpPositionAtLC_y = -1.*it_lc.z()*TMath::Tan(2.*(TMath::ATan(TMath::Exp(cp_eta))))*TMath::Sin(cp_phi);

        float xdist = cpPositionAtLC_x - recHitTools_.getPosition(detid_).x();
        float ydist = cpPositionAtLC_y - recHitTools_.getPosition(detid_).y();

        dist_eta_cp_lc->Fill(cp_eta - it_lc.eta());
        dist_phi_cp_lc->Fill(cp_phi - it_lc.phi());

        dist_x_cp_lc->Fill(cpPositionAtLC_x - it_lc.x());
        dist_y_cp_lc->Fill(cpPositionAtLC_y - it_lc.y());

        /*

        std::cout << "Layer: " << layer_ <<std::endl;     

        float edist = cp_eta - it_lc.eta();
        float phidist = cp_phi - it_lc.phi();
        float drdist = getDr(it_lc.eta(),it_lc.phi(),cp_eta,cp_phi);

        std::cout<<"Eta distance: " << edist <<endl;
        std::cout<<"phi distance: " << phidist <<endl;
        std::cout<<"dr distance: " << drdist <<endl;

        
        dist_dr_sim_lc->Fill(drdist);
        dist_eta_sim_lc->Fill(edist);
        dist_phi_sim_lc->Fill(phidist);
        dist_eta_phi_sim_lc->Fill(edist,phidist);
      

        */

       // if (tmprechits_[i0] != hitMap.end()){
        if(recHitTools_.isSilicon(tmprechits_[idx_[i0]])){
          hit_layer_lc_si[layer_-1]++;

          //dist_dr_cp_lc_sc->Fill(drdist);
          dist_eta_cp_lc_si->Fill(edist);
          dist_phi_cp_lc_si->Fill(phidist);
          dist_eta_phi_cp_lc_si->Fill(edist,phidist); 

          dist_x_cp_lc_si->Fill(xdist);
          dist_y_cp_lc_si->Fill(ydist);
          dist_x_y_cp_lc_si->Fill(xdist,ydist);
        }
        if(recHitTools_.isScintillator(tmprechits_[idx_[i0]])){
          hit_layer_lc_sc[layer_-1]++;

          //dist_dr_cp_lc_sc->Fill(drdist);
          dist_eta_cp_lc_sc->Fill(edist);
          dist_phi_cp_lc_sc->Fill(phidist);
          dist_eta_phi_cp_lc_sc->Fill(edist,phidist); 

          dist_x_cp_lc_sc->Fill(xdist);
          dist_y_cp_lc_sc->Fill(ydist);
          dist_x_y_cp_lc_sc->Fill(xdist,ydist);
        }
        hit_layer_lc[layer_-1]++;
       // }
      }
    }
  }
  
  


  cout<<"nTrack = "<<nTrack<<endl;
  //cout<<"nCaloparticles = "<<nCalo<<endl;
  //cout<<"nSimparticles = "<<nSim<<endl;
  //cout<<"nLCparticles = "<<nLC<<endl;
  cout<<"nSimhits = "<<nSimhits<<endl;
  cout<<"nRechits = "<<nRechits<<endl;
  //cout << "Num of LCs = "<< lcs.size() <<endl ;
  //cout<<"nLChits = "<<nLChits<<endl;


  hit_sim = nSimhits;
  hit_rec = nRechits;
  hit_lc = nLChits;

  hit_sim_histo->Fill(nSimhits);
  hit_rec_histo->Fill(nRechits);
  hit_lc_histo->Fill(nLChits);

  val_particle__sim_histo->Fill(nSim);
  val_particle__lc_histo->Fill(nLC);

  int lc_connectivity_counter = 0;
  int reccluster_connectivity_counter = 0;
  int simcluster_connectivity_counter = 0;

  int lc_si_connectivity_counter = 0;
  int reccluster_si_connectivity_counter = 0;
  int simcluster_si_connectivity_counter = 0;

  int lc_sc_connectivity_counter = 0;
  int reccluster_sc_connectivity_counter = 0;
  int simcluster_sc_connectivity_counter = 0;

  int max_lc_connectivity_counter = 0;
  int max_reccluster_connectivity_counter = 0;
  int max_simcluster_connectivity_counter = 0;

  int max_lc_si_connectivity_counter = 0;
  int max_reccluster_si_connectivity_counter = 0;
  int max_simcluster_si_connectivity_counter = 0;

  int max_lc_sc_connectivity_counter = 0;
  int max_reccluster_sc_connectivity_counter = 0;
  int max_simcluster_sc_connectivity_counter = 0;

  int miss_lc_connectivity_counter = 0;
  int miss_reccluster_connectivity_counter = 0;
  int miss_simcluster_connectivity_counter = 0;

  int miss_lc_si_connectivity_counter = 0;
  int miss_reccluster_si_connectivity_counter = 0;
  int miss_simcluster_si_connectivity_counter = 0;

  int miss_lc_sc_connectivity_counter = 0;
  int miss_reccluster_sc_connectivity_counter = 0;
  int miss_simcluster_sc_connectivity_counter = 0;

  int skip_lc_connectivity_counter = 0;
  int skip_reccluster_connectivity_counter = 0;
  int skip_simcluster_connectivity_counter = 0;

  int skip_lc_si_connectivity_counter = 0;
  int skip_reccluster_si_connectivity_counter = 0;
  int skip_simcluster_si_connectivity_counter = 0;

  int skip_lc_sc_connectivity_counter = 0;
  int skip_reccluster_sc_connectivity_counter = 0;
  int skip_simcluster_sc_connectivity_counter = 0;

  std::vector<int> skip_counters(9,0);


  std::vector<int> connectivity_counters = {lc_connectivity_counter, reccluster_connectivity_counter, simcluster_connectivity_counter,
                               lc_si_connectivity_counter, reccluster_si_connectivity_counter, simcluster_si_connectivity_counter,
                               lc_sc_connectivity_counter, reccluster_sc_connectivity_counter, simcluster_sc_connectivity_counter};

  std::vector<int> connectivity_miss_counters = {miss_lc_connectivity_counter, miss_reccluster_connectivity_counter, miss_simcluster_connectivity_counter,
                                   miss_lc_si_connectivity_counter, miss_reccluster_si_connectivity_counter, miss_simcluster_si_connectivity_counter,
                                   miss_lc_sc_connectivity_counter, miss_reccluster_sc_connectivity_counter, miss_simcluster_sc_connectivity_counter};

  std::vector<int> connectivity_max_counters = {max_lc_connectivity_counter, max_reccluster_connectivity_counter, max_simcluster_connectivity_counter,
                                   max_lc_si_connectivity_counter, max_reccluster_si_connectivity_counter, max_simcluster_si_connectivity_counter,
                                   max_lc_sc_connectivity_counter, max_reccluster_sc_connectivity_counter, max_simcluster_sc_connectivity_counter};

  std::vector<int> connectivity_skip_counters = {skip_lc_connectivity_counter, skip_reccluster_connectivity_counter, skip_simcluster_connectivity_counter,
                               skip_lc_si_connectivity_counter, skip_reccluster_si_connectivity_counter, skip_simcluster_si_connectivity_counter,
                               skip_lc_sc_connectivity_counter, skip_reccluster_sc_connectivity_counter, skip_simcluster_sc_connectivity_counter};


  std::vector<TH1F*> connectivity_histos = {connectivity_lc_histo, connectivity_rec_histo, connectivity_sim_histo,
                               connectivity_lc_si_histo, connectivity_rec_si_histo, connectivity_sim_si_histo,
                               connectivity_lc_sc_histo, connectivity_rec_sc_histo, connectivity_sim_sc_histo};

  std::vector<TH1F*> connectivity_w_histos = {connectivity_w_lc_histo, connectivity_w_rec_histo, connectivity_w_sim_histo,
                               connectivity_w_lc_si_histo, connectivity_w_rec_si_histo, connectivity_w_sim_si_histo,
                               connectivity_w_lc_sc_histo, connectivity_w_rec_sc_histo, connectivity_w_sim_sc_histo};

  std::vector<TH1F*> connectivity_max_histos = {connectivity_max_lc_histo, connectivity_max_rec_histo, connectivity_max_sim_histo,
                               connectivity_max_lc_si_histo, connectivity_max_rec_si_histo, connectivity_max_sim_si_histo,
                               connectivity_max_lc_sc_histo, connectivity_max_rec_sc_histo, connectivity_max_sim_sc_histo};

  std::vector<TH1F*> connectivity_miss_histos = {connectivity_miss_lc_histo, connectivity_miss_rec_histo, connectivity_miss_sim_histo,
                               connectivity_miss_lc_si_histo, connectivity_miss_rec_si_histo, connectivity_miss_sim_si_histo,
                               connectivity_miss_lc_sc_histo, connectivity_miss_rec_sc_histo, connectivity_miss_sim_sc_histo};

  std::vector<TH1F*> connectivity_skip_histos = {connectivity_skip_lc_histo, connectivity_skip_rec_histo, connectivity_skip_sim_histo,
                               connectivity_skip_lc_si_histo, connectivity_skip_rec_si_histo, connectivity_skip_sim_si_histo,
                               connectivity_skip_lc_sc_histo, connectivity_skip_rec_sc_histo, connectivity_skip_sim_sc_histo};


  std::vector<TH1F*> det_bool_histos = {det_bool_lc_histo, det_bool_rec_histo, det_bool_sim_histo,
                               det_bool_lc_si_histo, det_bool_rec_si_histo, det_bool_sim_si_histo,
                               det_bool_lc_sc_histo, det_bool_rec_sc_histo, det_bool_sim_sc_histo};

  std::vector<std::vector<int>> hit_layer_counters = {hit_layer_lc, hit_layer_rec, hit_layer_sim,
                               hit_layer_lc_si, hit_layer_rec_si, hit_layer_sim_si,
                               hit_layer_lc_sc, hit_layer_rec_sc, hit_layer_sim_sc};


  for(int i=0; i<47;i++){
    
    hit_layer_sim_histo[i]->Fill(hit_layer_sim[i]);
    hit_layer_rec_histo[i]->Fill(hit_layer_rec[i]);
    hit_layer_lc_histo[i]->Fill(hit_layer_lc[i]);  

    hit_layer_sim_si_histo[i]->Fill(hit_layer_sim_si[i]);
    hit_layer_rec_si_histo[i]->Fill(hit_layer_rec_si[i]);
    hit_layer_lc_si_histo[i]->Fill(hit_layer_lc_si[i]);

    hit_layer_sim_sc_histo[i]->Fill(hit_layer_sim_sc[i]);
    hit_layer_rec_sc_histo[i]->Fill(hit_layer_rec_sc[i]);
    hit_layer_lc_sc_histo[i]->Fill(hit_layer_lc_sc[i]);  
    
    det_sim_histo->Fill(i,hit_layer_sim[i]);
    det_sim_si_histo->Fill(i,hit_layer_sim_si[i]);
    det_sim_sc_histo->Fill(i,hit_layer_sim_sc[i]);

    det_rec_histo->Fill(i,hit_layer_rec[i]);
    det_rec_si_histo->Fill(i,hit_layer_rec_si[i]);
    det_rec_sc_histo->Fill(i,hit_layer_rec_sc[i]);

    det_lc_histo->Fill(i,hit_layer_lc[i]);
    det_lc_si_histo->Fill(i,hit_layer_lc_si[i]);
    det_lc_sc_histo->Fill(i,hit_layer_lc_sc[i]);

    // Loop through the different counters for the various connectivity scores and fill the histograms

    for(unsigned j=0; j<connectivity_counters.size();j++){

      if(hit_layer_counters[j][i]!=0){
        det_bool_histos[j]->Fill(i);
        connectivity_counters[j]++;
        connectivity_skip_counters[j]++;
	skip_counters[j]=skip_;
        if(connectivity_miss_counters[j]!=0){connectivity_miss_histos[j]->Fill(connectivity_miss_counters[j]);}
        connectivity_miss_counters[j]=0;
        if(i==46){
          connectivity_histos[j]->Fill(connectivity_counters[j]);
          connectivity_skip_histos[j]->Fill(connectivity_skip_counters[j]);

          if(connectivity_max_counters[j]<connectivity_counters[j]){connectivity_max_counters[j] = connectivity_counters[j];}
          connectivity_max_histos[j]->Fill(connectivity_max_counters[j]);

          for(int k=i+1-connectivity_counters[j];k<i+1;k++){
            connectivity_w_histos[j]->Fill(k,connectivity_counters[j]);
          }
        }
      }

      else{

        connectivity_miss_counters[j]++;

        if(skip_counters[j] !=0){
          skip_counters[j]--;
          connectivity_skip_counters[j]++;
        }
        else{
          if(connectivity_skip_counters[j]-skip_>0){connectivity_skip_histos[j]->Fill(connectivity_skip_counters[j]-skip_);}

          connectivity_skip_counters[j] = 0;
        }
        if(connectivity_counters[j]!=0){connectivity_histos[j]->Fill(connectivity_counters[j]);}
        for(int k=i-connectivity_counters[j]; k<i;k++){
          connectivity_w_histos[j]->Fill(k,connectivity_counters[j]);
        }
        if(connectivity_max_counters[j]<connectivity_counters[j]){connectivity_max_counters[j] = connectivity_counters[j];}
        connectivity_counters[j]=0;
        if(i==46){
          connectivity_max_histos[j]->Fill(connectivity_max_counters[j]);
          connectivity_miss_histos[j] ->Fill(connectivity_miss_counters[j]);
        }
      }
    }


    /*

    // Create Plots for Simhits 

    // Silicon

    if (hit_layer_sim_si[i]!=0){
      det_bool_sim_si_histo->Fill(i);
      simcluster_si_connectivity_counter++;
      if (i==46){
        connectivity_sim_si_histo->Fill(simcluster_si_connectivity_counter);
        if (max_simcluster_si_connectivity_counter < simcluster_si_connectivity_counter){max_simcluster_si_connectivity_counter = simcluster_si_connectivity_counter;}
        connectivity_max_sim_si_histo->Fill(max_simcluster_si_connectivity_counter);
        for(int j=i+1-simcluster_si_connectivity_counter; j<i+1;j++){
          connectivity_w_sim_si_histo->Fill(j, simcluster_si_connectivity_counter);
        }
      }
    }
    else {
      connectivity_sim_si_histo->Fill(simcluster_si_connectivity_counter);
      for(int j=i-simcluster_si_connectivity_counter; j<i;j++){
        connectivity_w_sim_si_histo->Fill(j, simcluster_si_connectivity_counter);
      }
      if (max_simcluster_si_connectivity_counter < simcluster_si_connectivity_counter){max_simcluster_si_connectivity_counter = simcluster_si_connectivity_counter;}
      simcluster_si_connectivity_counter = 0;
    }

    // Scintillator

    if (hit_layer_sim_sc[i]!=0){
      det_bool_sim_sc_histo->Fill(i);
      simcluster_sc_connectivity_counter++;
      if (i==46){
        connectivity_sim_sc_histo->Fill(simcluster_sc_connectivity_counter);
        if (max_simcluster_sc_connectivity_counter < simcluster_sc_connectivity_counter){max_simcluster_sc_connectivity_counter = simcluster_sc_connectivity_counter;}
        connectivity_max_sim_sc_histo->Fill(max_simcluster_sc_connectivity_counter);
        for(int j=i+1-simcluster_sc_connectivity_counter; j<i+1;j++){
          connectivity_w_sim_sc_histo->Fill(j, simcluster_sc_connectivity_counter);
        }
      }
    }
    else {
      connectivity_sim_sc_histo->Fill(simcluster_sc_connectivity_counter);
      for(int j=i-simcluster_sc_connectivity_counter; j<i;j++){
        connectivity_w_sim_sc_histo->Fill(j, simcluster_sc_connectivity_counter);
      }
      if (max_simcluster_sc_connectivity_counter < simcluster_sc_connectivity_counter){max_simcluster_sc_connectivity_counter = simcluster_sc_connectivity_counter;}
      simcluster_sc_connectivity_counter = 0;
    }

    // Total

    if (hit_layer_sim[i]!=0){
      det_bool_sim_histo->Fill(i);
      simcluster_connectivity_counter++;
      if (i==46){
        connectivity_sim_histo->Fill(simcluster_connectivity_counter);
        if (max_simcluster_connectivity_counter < simcluster_connectivity_counter){max_simcluster_connectivity_counter = simcluster_connectivity_counter;}
        connectivity_max_sim_histo->Fill(max_simcluster_connectivity_counter);
        for(int j=i+1-simcluster_connectivity_counter; j<i+1;j++){
          connectivity_w_sim_histo->Fill(j, simcluster_connectivity_counter);
        }
      }
    }
    else {
      connectivity_sim_histo->Fill(simcluster_connectivity_counter);
      for(int j=i-simcluster_connectivity_counter; j<i;j++){
        connectivity_w_sim_histo->Fill(j, simcluster_connectivity_counter);
      }
      if (max_simcluster_connectivity_counter < simcluster_connectivity_counter){max_simcluster_connectivity_counter = simcluster_connectivity_counter;}
      simcluster_connectivity_counter = 0;
    }

    // Create Plots for Rechits

    // Si

    if (hit_layer_rec_si[i]!=0){
      det_bool_rec_si_histo->Fill(i);
      reccluster_si_connectivity_counter++;
      if (i==46){
        connectivity_rec_si_histo->Fill(reccluster_si_connectivity_counter);
        if (max_reccluster_si_connectivity_counter < reccluster_si_connectivity_counter){max_reccluster_si_connectivity_counter = reccluster_si_connectivity_counter;}
        connectivity_max_rec_si_histo->Fill(max_reccluster_si_connectivity_counter);
        for(int j=i+1-reccluster_si_connectivity_counter; j<i+1;j++){
          connectivity_w_rec_si_histo->Fill(j, reccluster_si_connectivity_counter);
        }
      }
    }
    else {
      connectivity_rec_si_histo->Fill(reccluster_si_connectivity_counter);
      for(int j=i-reccluster_si_connectivity_counter; j<i;j++){
        connectivity_w_rec_si_histo->Fill(j, reccluster_si_connectivity_counter);
      }
      if (max_reccluster_si_connectivity_counter < reccluster_si_connectivity_counter){max_reccluster_si_connectivity_counter = reccluster_si_connectivity_counter;}
      reccluster_si_connectivity_counter = 0;
    }

    // Scintillator

    if (hit_layer_rec_sc[i]!=0){
      det_bool_rec_sc_histo->Fill(i);
      reccluster_sc_connectivity_counter++;
      if (i==46){
        connectivity_rec_sc_histo->Fill(reccluster_sc_connectivity_counter);
        if (max_reccluster_sc_connectivity_counter < reccluster_sc_connectivity_counter){max_reccluster_sc_connectivity_counter = reccluster_sc_connectivity_counter;}
        connectivity_max_rec_sc_histo->Fill(max_reccluster_sc_connectivity_counter);
        for(int j=i+1-reccluster_sc_connectivity_counter; j<i+1;j++){
          connectivity_w_rec_sc_histo->Fill(j, reccluster_sc_connectivity_counter);
        }
      }
    }
    else {
      connectivity_rec_sc_histo->Fill(reccluster_sc_connectivity_counter);
      for(int j=i-reccluster_sc_connectivity_counter; j<i;j++){
        connectivity_w_rec_sc_histo->Fill(j, reccluster_sc_connectivity_counter);
      }
      if (max_reccluster_sc_connectivity_counter < reccluster_sc_connectivity_counter){max_reccluster_sc_connectivity_counter = reccluster_sc_connectivity_counter;}
      reccluster_sc_connectivity_counter = 0;
    }
    
    // Total

    if (hit_layer_rec[i]!=0){
      det_bool_rec_histo->Fill(i);
      reccluster_connectivity_counter++;
      if (i==46){
        connectivity_rec_histo->Fill(reccluster_connectivity_counter);
        if (max_reccluster_connectivity_counter < reccluster_connectivity_counter){max_reccluster_connectivity_counter = reccluster_connectivity_counter;}
        connectivity_max_rec_histo->Fill(max_reccluster_connectivity_counter);
        for(int j=i+1-reccluster_connectivity_counter; j<i+1;j++){
          connectivity_w_rec_histo->Fill(j, reccluster_connectivity_counter);
        }
      }
    }
    else {
      connectivity_rec_histo->Fill(reccluster_connectivity_counter);
      for(int j=i-reccluster_connectivity_counter; j<i;j++){
        connectivity_w_rec_histo->Fill(j, reccluster_connectivity_counter);
      }
      if (max_reccluster_connectivity_counter < reccluster_connectivity_counter){max_reccluster_connectivity_counter = reccluster_connectivity_counter;}
      reccluster_connectivity_counter = 0;
    }

    // Create Plots for LC

    // Si

    if (hit_layer_lc_si[i]!=0){
      det_bool_lc_si_histo->Fill(i);
      lc_si_connectivity_counter++;
      if (i==46){
        connectivity_lc_si_histo->Fill(lc_si_connectivity_counter);
        if (max_lc_si_connectivity_counter < lc_si_connectivity_counter){max_lc_si_connectivity_counter = lc_si_connectivity_counter;}
        connectivity_max_lc_si_histo->Fill(max_lc_si_connectivity_counter);
        for(int j=i+1-lc_si_connectivity_counter; j<i+1;j++){
          connectivity_w_lc_si_histo->Fill(j, lc_si_connectivity_counter);
          }
        }
      }
    else {
      connectivity_lc_si_histo->Fill(lc_si_connectivity_counter);
      for(int j=i-lc_si_connectivity_counter; j<i;j++){
        connectivity_w_lc_si_histo->Fill(j, lc_si_connectivity_counter);
        }
      if (max_lc_si_connectivity_counter < lc_si_connectivity_counter){max_lc_si_connectivity_counter = lc_si_connectivity_counter;}
      lc_si_connectivity_counter = 0;
      }

    // Scintillator

    if (hit_layer_lc_sc[i]!=0){
      det_bool_lc_sc_histo->Fill(i);
      lc_sc_connectivity_counter++;
      if (i==46){
        connectivity_lc_sc_histo->Fill(lc_sc_connectivity_counter);
        if (max_lc_sc_connectivity_counter < simcluster_sc_connectivity_counter){max_lc_sc_connectivity_counter = lc_sc_connectivity_counter;}
        connectivity_max_lc_sc_histo->Fill(max_lc_sc_connectivity_counter);
        for(int j=i+1-lc_sc_connectivity_counter; j<i+1;j++){
          connectivity_w_lc_sc_histo->Fill(j, lc_sc_connectivity_counter);
          }
        }
      }
    else {
      connectivity_lc_sc_histo->Fill(lc_sc_connectivity_counter);
      for(int j=i-lc_sc_connectivity_counter; j<i;j++){
        connectivity_w_lc_sc_histo->Fill(j, lc_sc_connectivity_counter);
        }
      if (max_lc_sc_connectivity_counter < lc_sc_connectivity_counter){max_lc_si_connectivity_counter = lc_sc_connectivity_counter;}
      lc_sc_connectivity_counter = 0;
      }


    if (hit_layer_lc[i]!=0){
      det_bool_lc_histo->Fill(i);
      lc_connectivity_counter++;
      if (i==46){
        connectivity_lc_histo->Fill(lc_connectivity_counter);
        if (max_lc_connectivity_counter < lc_connectivity_counter){max_lc_connectivity_counter = lc_connectivity_counter;}
        connectivity_max_lc_histo->Fill(max_lc_connectivity_counter);
        for(int j=i+1-lc_connectivity_counter; j<i+1;j++){
          connectivity_w_lc_histo->Fill(j, lc_connectivity_counter);
          }
        }
      }
    else {
      connectivity_lc_histo->Fill(lc_connectivity_counter);
      for(int j=i-lc_connectivity_counter; j<i;j++){
        connectivity_w_lc_histo->Fill(j, lc_connectivity_counter);
        }
      if (max_lc_connectivity_counter < lc_connectivity_counter){max_lc_connectivity_counter = lc_connectivity_counter;}
      lc_connectivity_counter = 0;
      }

    */
    }

  tree->Fill();


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
