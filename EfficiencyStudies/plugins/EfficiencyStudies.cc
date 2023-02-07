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
#include <any>
#include <iomanip>
#include <cmath>
#include <typeinfo>

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

  std::vector<std::string> detectors;
  std::vector<std::string> objects;
  hgcal::RecHitTools recHitTools_;

  // ----------member data ---------------------------
  edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
  edm::EDGetTokenT<std::vector<CaloParticle>> caloParticlesToken_;
  edm::EDGetTokenT<std::vector<ticl::Trackster>> ticlTrackstersMergeToken;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsEEToken_; 
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsFHToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsBHToken_;
  edm::EDGetTokenT<reco::CaloClusterCollection> hgcalLayerClustersToken_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;

  TTree *tree = new TTree("tree","tree");

  typedef std::map<std::string, std::vector<TH1F*>> detmap;
  typedef std::map<std::string, detmap> objmap;

  std::map<std::string, objmap> plot_map;

  typedef std::map<std::string, std::vector<TH2F*>> detmap2D;
  typedef std::map<std::string, detmap2D> objmap2D;

  std::map<std::string, objmap2D> plot_2D_map;

  // double trackPtMin_; 
  int skip_;
  std::string eta_;

  // Energy

  TH1F *energy_rechits_histo;

  // Missing Simhits

  TH2F *missing_simhits_histo;

  // Validate Particles

  TH1F *val_particle__sim_histo;
  TH1F *val_particle__rec_histo;
  TH1F *val_particle__lc_histo;

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
      ticlTrackstersMergeToken(consumes<std::vector<ticl::Trackster> >(iConfig.getParameter<edm::InputTag>("Tracksters"))),
      hgcalRecHitsEEToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsEE"))),
      hgcalRecHitsFHToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsFH"))),
      hgcalRecHitsBHToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsBH"))),
      hgcalLayerClustersToken_(consumes<reco::CaloClusterCollection>(iConfig.getParameter<edm::InputTag>("hgcalLayerClusters"))),
      // trackPtMin_(iConfig.getParameter<double>("trackPtMin")), 
      caloGeomToken_(esConsumes<CaloGeometry, CaloGeometryRecord>()),
      skip_(iConfig.getParameter<int>("skip")),
      eta_(iConfig.getParameter<std::string>("eta")){

  detectors = {"", "Si", "Si 120", "Si 200", "Si 300", "Sc"};
  objects = {"Simhits", "Rechits", "LCs"};



#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif

  //recHitTools_.reset(new hgcal::RecHitTools());
  //now do what ever initialization is needed

  usesResource("TFileService");
  edm::Service<TFileService> file;

  tree = file->make<TTree>("tree","hgc analyzer");

  // One Off Histograms

  // Energy histograms

  plot_map["Energy"]["Rechits"][""].push_back(new TH1F("Energy Rechits", "Energy Rechits", 100,0,5));

  // Missing Hits 2D Histogram

  missing_simhits_histo = new TH2F("Missing Simhits","Missing Simhits",47,0,47,500,1,501);

  // Extensive Histograms

  TString t_name;
  for(const auto& obj : objects){

    // Hit Plots

    t_name = obj;
    plot_map["Hit Plots"][obj][""].push_back(new TH1F(t_name, t_name, 500, 0, 500));

    // Validation Plots

    t_name = obj + " Cluster";
    plot_map["Validation Plots"][obj][""].push_back(new TH1F(t_name, t_name, 500, 0, 500));

    for(const auto& det : detectors){

      // Hit Plots

      t_name = "Detector Hit Plot: " + obj + " " + det;     
      plot_map["Detector Hit Plots"][obj][det].push_back(new TH1F(t_name, t_name, 47, 0, 47));

      t_name = "Detector Bool Hit Plot: " + obj + " " + det; 
      plot_map["Detector Bool Hit Plots"][obj][det].push_back(new TH1F(t_name, t_name, 47,0,47));

      // Connectivity Plots
    
      t_name = "Connectivity Plot: " + obj + " " + det;
      plot_map["Connectivity Plots"][obj][det].push_back(new TH1F(t_name, t_name, 48,0,48));

      t_name = "W Connectivity Plot: " + obj + " " + det;
      plot_map["W Connectivity Plots"][obj][det].push_back(new TH1F(t_name, t_name, 48,0,48));

      t_name = "Skip Connectivity Plot: " + obj + " " + det;
      plot_map["Skip Connectivity Plots"][obj][det].push_back(new TH1F(t_name, t_name, 48,0,48));

      t_name = "Max Connectivity Plot: " + obj + " " + det;
      plot_map["Max Connectivity Plots"][obj][det].push_back(new TH1F(t_name, t_name, 48,0,48));

      t_name = "Miss Connectivity Plot: " + obj + " " + det;
      plot_map["Miss Connectivity Plots"][obj][det].push_back(new TH1F(t_name, t_name, 48,0,48));

      // Diff Plots

      t_name = "Diff Eta Plots: "  + obj + " " + det;
      plot_map["Diff Eta Plots"][obj][det].push_back(new TH1F(t_name, t_name, 300,-0.15,0.15));

      t_name = "Diff Phi Plots: "  + obj + " " + det;
      plot_map["Diff Phi Plots"][obj][det].push_back(new TH1F(t_name, t_name, 300,-0.15,0.15));

      t_name = "Diff x Plots: "  + obj + " " + det;
      plot_map["Diff x Plots"][obj][det].push_back(new TH1F(t_name, t_name, 200,-10,10));

      t_name = "Diff y Plots: "  + obj + " " + det;
      plot_map["Diff y Plots"][obj][det].push_back(new TH1F(t_name, t_name, 200,-10,10));

      // 2D Diff Plots

      t_name = "Diff Eta Phi Plots: "  + obj + " " + det;
      plot_2D_map["Diff Eta Phi Plots"][obj][det].push_back(new TH2F(t_name,t_name,200,-0.2,0.2,200,-0.2,0.2));

      t_name = "Diff x y Plots: "  + obj + " " + det;
      plot_2D_map["Diff x y Plots"][obj][det].push_back(new TH2F(t_name,t_name,200,-10,10,200,-10,10));

      for(int i=0; i<47; i++){
        std::stringstream nlayer;
        nlayer << i;

        t_name = "Hit Layer Plot: " + obj + " " + det + " Layer_" + nlayer.str();    
        plot_map["Hit Layer Plots"][obj][det].push_back(new TH1F(t_name,t_name,20,0,20));

        t_name = "Diff Eta Layer Plot: " + obj + " " + det + " Layer_" + nlayer.str(); 
        plot_map["Diff Eta Layer Plots"][obj][det].push_back(new TH1F(t_name,t_name,200,-10,10));

        t_name = "Diff Phi Layer Plot: " + obj + " " + det + " Layer_" + nlayer.str();
        plot_map["Diff Phi Layer Plots"][obj][det].push_back(new TH1F(t_name, t_name,200,-10,10));

        t_name = "Diff x Layer Plot: " + obj + " " + det + " Layer_" + nlayer.str();
        plot_map["Diff x Layer Plots"][obj][det].push_back(new TH1F(t_name,t_name,200,-0,0.2));

        t_name = "Diff y Layer Plot: " + obj + " " + det + " Layer_" + nlayer.str();       
        plot_map["Diff y Layer Plots"][obj][det].push_back(new TH1F(t_name,t_name,200,-0.2,0.2));
      }
    }
  }

}

EfficiencyStudies::~EfficiencyStudies() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty


  //TLegend *leg = new TLegend(0.68,0.72,0.98,0.92); // Is dynamic memory allocation necessary?

  std::cout << "Write stuff" <<std::endl;

  plot_map["Energy"]["Rechits"][""].front()->Write();
  missing_simhits_histo->Write();
  for(auto& obj: objects){
    plot_map["Hit Plots"][obj][""].front()->Write();
    plot_map["Validation Plots"][obj][""].front()->Write();
      
    for(auto& det: detectors){
      plot_map["Detector Hit Plots"][obj][det].front()->Write();
      plot_map["Detector Bool Hit Plots"][obj][det].front()->Write();
      plot_map["Connectivity Plots"][obj][det].front()->Write();
      plot_map["W Connectivity Plots"][obj][det].front()->Write();
      plot_map["Skip Connectivity Plots"][obj][det].front()->Write();
      plot_map["Max Connectivity Plots"][obj][det].front()->Write();
      plot_map["Miss Connectivity Plots"][obj][det].front()->Write();
      plot_map["Diff Eta Plots"][obj][det].front()->Write();
      plot_map["Diff Phi Plots"][obj][det].front()->Write();
      plot_map["Diff x Plots"][obj][det].front()->Write();
      plot_map["Diff y Plots"][obj][det].front()->Write();
      plot_2D_map["Diff Eta Phi Plots"][obj][det].front()->Write();
      plot_2D_map["Diff x y Plots"][obj][det].front()->Write();
    }
  }

  /*
  TCanvas *c1 = new TCanvas("c1","c1"); // Is dynamic memory allocation necessary?

  std::vector<int> colors = {1, 4, 2};
  TString folder = "/eos/user/m/mmatthew/www/EfficiencyStudies/Analyzer/"+eta_+"/";
  std::vector<float> pos = {0.70,0.75,0.9,0.9};
  std::vector<TString> axes;

  // Energy Histograms
  std::cout<<"Energy"<<std::endl;
  axes={"Occurence", "Energy"};

  createTH1Plot(c1, plot_map["Energy"]["Rechits"][""].front(), "EnergySimhits.png", axes, folder);

  // Missing Hits Plots
  std::cout<<"Missing Hits"<<std::endl;
  axes = {"Layer", "Events"};
  createTH2Plot(c1, missing_simhits_histo, "MissingSimhits.png", axes, folder);

  // Validation Plots

  std::cout<<"Validation"<<std::endl;
  axes = {"# Particles", "# Events"};
  createTH1Plot(c1, plot_map["Validation Plots"]["Simhits"][""].front(), "SimParticles.png", axes, folder);

  axes = {"# Clusters","# Occurence"};
  createTH1Plot(c1, plot_map["Validation Plots"]["LCs"][""].front(), "LayerClusters.png", axes, folder);

  // Hit Plots

  axes = {"Hits", "Events"};

  std::cout<<"Hit"<<std::endl;
  createTH1Plot(c1,plot_map["Hit Plots"]["Simhits"][""].front(), "Simhits.png",axes,folder);
  createTH1Plot(c1,plot_map["Hit Plots"]["Rechits"][""].front(), "Rechits.png",axes,folder);
  createTH1Plot(c1,plot_map["Hit Plots"]["LCs"][""].front(), "LCs.png",axes,folder);

  std::vector<std::string> objects = {"Simhits", "Rechits", "LCs"};
  std::vector<std::string> detectors = {"", "Si", "Si 120", "Si 200", "Si 300", "Sc"};
  std::vector<TString> legend =  {"SimHits", "RecHits", "LCs"};
  std::vector<TH1F*> hists = {};
  TString t_name;
  TString ttitle;

  for(const auto& det : detectors){


    // Det Hits Plots


    axes = {"Layer","# Clusters"};
    std::cout<<"Detector Bool Hit Plots"<<std::endl;
    createTH1Plot(c1, plot_map["Detector Hit Plots"]["Simhits"][det].front(),"Layerwise_Simhits_" + det + ".png", axes, folder);
    createTH1Plot(c1, plot_map["Detector Hit Plots"]["Rechits"][det].front(),"Layerwise_Rechits_" + det + ".png", axes, folder);
    createTH1Plot(c1, plot_map["Detector Hit Plots"]["LCs"][det].front(),"Layerwise_LCs_" + det + ".png", axes, folder);


    // Connectivity Plots

    axes = {"Chain Length", "Occurence"};
    c1->SetLogy();

    createTH1Plot(c1, plot_map["Connectivity Plots"]["Simhits"][det].front(),"Connectivity_Simhits_" + det + ".png", axes, folder);
    createTH1Plot(c1, plot_map["Connectivity Plots"]["Rechits"][det].front(),"Connectivity_Rechits_" + det + ".png", axes, folder);
    createTH1Plot(c1, plot_map["Connectivity Plots"]["LCs"][det].front(),"Connectivity_LC_" + det + ".png", axes, folder);

    plot_map["Connectivity Plots"]["Simhits"][det].front() ->SetStats(0);
    plot_map["Connectivity Plots"]["Rechits"][det].front() ->SetStats(0);
    plot_map["Connectivity Plots"]["LCs"][det].front() ->SetStats(0);

    t_name = "Connectivity_Stack_" + det + ".png";
    hists = {plot_map["Connectivity Plots"]["Simhits"][det].front(), plot_map["Connectivity Plots"]["Simhits"][det].front(), plot_map["Connectivity Plots"]["Simhits"][det].front()};
    createTHSPlot(c1, colors, hists, t_name, axes, legend, {0.60,0.73,0.8,0.88}, folder, "Total Connectivity: " + det); //{0.60,0.70,0.9,0.9}

    // Skip Connectivity Plots

    std::stringstream ss;
    ss << "Connectivity_Skip_" <<std::setw(1) << std::setfill('1') << skip_ << "_" << det;
    t_name = ss.str();

    createTH1Plot(c1, plot_map["Skip Connectivity Plots"]["Simhits"][det].front(),t_name+"Simhits_" + det + ".png", axes, folder);
    createTH1Plot(c1, plot_map["Skip Connectivity Plots"]["Rechits"][det].front(),t_name+"Rechits_" + det + ".png", axes, folder);
    createTH1Plot(c1, plot_map["Skip Connectivity Plots"]["LCs"][det].front(),t_name+"LC_" + det + ".png", axes, folder);

    plot_map["Skip Connectivity Plots"]["Simhits"][det].front() ->SetStats(0);
    plot_map["Skip Connectivity Plots"]["Rechits"][det].front() ->SetStats(0);
    plot_map["Skip Connectivity Plots"]["LCs"][det].front() ->SetStats(0);

    hists = {plot_map["Skip Connectivity Plots"]["Simhits"][det].front(), plot_map["Skip Connectivity Plots"]["Rechits"][det].front(), plot_map["Skip Connectivity Plots"]["LCs"][det].front()};
    createTHSPlot(c1, colors, hists, t_name + "Stack.png", axes, legend, {0.60,0.73,0.8,0.88}, folder, "Skip Connectivity: " + det);

    //Max

    c1->SetLogy(0);

    createTH1Plot(c1, plot_map["Max Connectivity Plots"]["Simhits"][det].front(),"Connectivity_Max_SimCluster_" + det +  ".png", axes, folder);
    createTH1Plot(c1, plot_map["Max Connectivity Plots"]["Rechits"][det].front(),"Connectivity_Max_RecCluster_" + det +  ".png", axes, folder);
    createTH1Plot(c1, plot_map["Max Connectivity Plots"]["LCs"][det].front(),"Connectivity_Max_LCs_" + det +  ".png", axes, folder);

    plot_map["Max Connectivity Plots"]["Simhits"][det].front() ->SetStats(0);
    plot_map["Max Connectivity Plots"]["Rechits"][det].front() ->SetStats(0);
    plot_map["Max Connectivity Plots"]["LCs"][det].front() ->SetStats(0);

    t_name = "Connectivity_Max_Stack_" + det + ".png";
    hists = {plot_map["Max Connectivity Plots"]["Simhits"][det].front(), plot_map["Max Connectivity Plots"]["Rechits"][det].front(), plot_map["Max Connectivity Plots"]["LCs"][det].front()};
    createTHSPlot(c1, colors, hists, t_name, axes, legend, {0.40,0.73,0.6,0.88}, folder, "Max Connectivity: " + det);

    // Miss

    createTH1Plot(c1, plot_map["Miss Connectivity Plots"]["Simhits"][det].front(),"Connectivity_Miss_SimCluster_" + det +  ".png", axes, folder);
    createTH1Plot(c1, plot_map["Miss Connectivity Plots"]["Rechits"][det].front(),"Connectivity_Miss_RecCluster_" + det +  ".png", axes, folder);
    createTH1Plot(c1, plot_map["Miss Connectivity Plots"]["LCs"][det].front(),"Connectivity_Miss_LCs_" + det +  ".png", axes, folder);

    plot_map["Miss Connectivity Plots"]["Simhits"][det].front() ->SetStats(0);
    plot_map["Miss Connectivity Plots"]["Rechits"][det].front() ->SetStats(0);
    plot_map["Miss Connectivity Plots"]["LCs"][det].front() ->SetStats(0);

    t_name = "Connectivity_Miss_Stack_" + det + ".png";
    hists = {plot_map["Miss Connectivity Plots"]["Simhits"][det].front(), plot_map["Miss Connectivity Plots"]["Rechits"][det].front(), plot_map["Miss Connectivity Plots"]["LCs"][det].front()};
    createTHSPlot(c1, colors, hists, t_name, axes, legend, pos, folder, "Miss Connectivity: " + det);

    // Weighted

    createTH1Plot(c1, plot_map["W Connectivity Plots"]["Simhits"][det].front(),"Connectivity_W_SimCluster_" + det +  ".png", axes, folder);
    createTH1Plot(c1, plot_map["W Connectivity Plots"]["Rechits"][det].front(),"Connectivity_W_RecCluster_" + det +  ".png", axes, folder);
    createTH1Plot(c1, plot_map["W Connectivity Plots"]["LCs"][det].front(),"Connectivity_W_LCs_" + det +  ".png", axes, folder);

    plot_map["W Connectivity Plots"]["Simhits"][det].front() ->SetStats(0);
    plot_map["W Connectivity Plots"]["Rechits"][det].front() ->SetStats(0);
    plot_map["W Connectivity Plots"]["LCs"][det].front() ->SetStats(0);

    t_name = "Connectivity_Weighted_Stack_" + det + ".png";
    hists = {plot_map["W Connectivity Plots"]["Simhits"][det].front(), plot_map["W Connectivity Plots"]["Rechits"][det].front(), plot_map["W Connectivity Plots"]["LCs"][det].front()};
    createTHSPlot(c1, colors, hists, t_name, axes, legend, pos, folder, "W Connectivity: " + det);

    // Efficiency Plots (all hits)

    axes = {"Layer", "Efficiency"};

    plot_map["Detector Hit Plots"]["Simhits"][det].front()->SetStats(0);
    TH1F *eff_sim_rec_histo = (TH1F*)plot_map["Detector Hit Plots"]["Rechits"][det].front()->Clone("Efficiency Sim vs Rec");
    eff_sim_rec_histo->SetLineColor(kBlue);
    eff_sim_rec_histo->SetStats(0);
    TH1F *eff_sim_lc_histo = (TH1F*)plot_map["Detector Hit Plots"]["LCs"][det].front()->Clone("Efficiency Sim vs LC");
    eff_sim_lc_histo->SetLineColor(kRed);
    eff_sim_lc_histo->SetStats(0);

    eff_sim_rec_histo->Divide(plot_map["Detector Hit Plots"]["Simhits"][det].front());
    eff_sim_lc_histo->Divide(plot_map["Detector Hit Plots"]["Simhits"][det].front());

    hists = {eff_sim_rec_histo, eff_sim_lc_histo};
    createTHSPlot(c1, {4,2}, hists, "Efficiency_Plots_"+det+".png", axes, {"RecHits", "LCs"}, {0.15,0.15,0.35,0.25}, folder, "Efficiency: " + det);\

    // Efficiency Plots (boolean)

    axes = {"Layer", "Efficiency"};

    plot_map["Detector Bool Hit Plots"]["Simhits"][det].front()->SetStats(0);
    TH1F *bool_eff_sim_rec_histo = (TH1F*)plot_map["Detector Bool Hit Plots"]["Rechits"][det].front()->Clone("Efficiency Sim vs Rec");
    bool_eff_sim_rec_histo->SetLineColor(kBlue);
    bool_eff_sim_rec_histo->SetStats(0);
    TH1F *bool_eff_sim_lc_histo = (TH1F*)plot_map["Detector Bool Hit Plots"]["LCs"][det].front()->Clone("Efficiency Sim vs LC");
    bool_eff_sim_lc_histo->SetLineColor(kRed);
    bool_eff_sim_lc_histo->SetStats(0);

    bool_eff_sim_rec_histo->Divide(plot_map["Detector Bool Hit Plots"]["Simhits"][det].front());
    bool_eff_sim_lc_histo->Divide(plot_map["Detector Bool Hit Plots"]["Simhits"][det].front());

    hists = {bool_eff_sim_rec_histo, bool_eff_sim_lc_histo};
    createTHSPlot(c1, {4,2}, hists, "Boolean_Efficiency_Plots_"+det+".png", axes, {"RecHits", "LCs"}, {0.15,0.15,0.35,0.25}, folder, "Boolean Efficiency: " + det);


    // Distance

    axes = {"#Delta #eta", "# Occurence"};
    createTH1Plot(c1, plot_map["Diff Eta Plots"]["Simhits"][det].front(),"Diff_Eta_CP_Sim_" + det + ".png", axes, folder);
    createTH1Plot(c1, plot_map["Diff Eta Plots"]["Rechits"][det].front(),"Diff_Eta_CP_Rec_" + det + ".png", axes, folder);
    createTH1Plot(c1, plot_map["Diff Eta Plots"]["LCs"][det].front(),"Diff_Eta_CP_LC_" + det + ".png", axes, folder);

    axes = {"#Delta #phi", "# Occurence"};
    createTH1Plot(c1, plot_map["Diff Phi Plots"]["Simhits"][det].front(),"Diff_Phi_CP_Sim_" + det + ".png", axes, folder);
    createTH1Plot(c1, plot_map["Diff Phi Plots"]["Rechits"][det].front(),"Diff_Phi_CP_Rec_" + det + ".png", axes, folder);
    createTH1Plot(c1, plot_map["Diff Phi Plots"]["LCs"][det].front(),"Diff_Phi_CP_LC_" + det + ".png", axes, folder);

    axes = {"#Delta x", "# Occurence"};
    createTH1Plot(c1, plot_map["Diff x Plots"]["Simhits"][det].front(),"Diff_x_CP_Sim_" + det + ".png", axes, folder);
    createTH1Plot(c1, plot_map["Diff x Plots"]["Rechits"][det].front(),"Diff_x_CP_Rec_" + det + ".png", axes, folder);
    createTH1Plot(c1, plot_map["Diff x Plots"]["LCs"][det].front(),"Diff_x_CP_LC_" + det + ".png", axes, folder);

    axes = {"#Delta y", "# Occurence"};
    createTH1Plot(c1, plot_map["Diff y Plots"]["Simhits"][det].front(),"Diff_Phi_CP_Sim_" + det + ".png", axes, folder);
    createTH1Plot(c1, plot_map["Diff y Plots"]["Rechits"][det].front(),"Diff_Phi_CP_Rec_" + det + ".png", axes, folder);
    createTH1Plot(c1, plot_map["Diff y Plots"]["LCs"][det].front(),"Diff_Phi_CP_LC_" + det + ".png", axes, folder);

    // 2D Diff Plots

    axes = {"#Delta #eta", "#Delta #phi"};
    createTH2Plot(c1, plot_2D_map["Diff Eta Phi Plots"]["Simhits"][det].front(),"Diff_Eta_Phi_CP_Sim_" + det + ".png", axes, folder);
    createTH2Plot(c1, plot_2D_map["Diff Eta Phi Plots"]["Rechits"][det].front(),"Diff_Eta_Phi_CP_Rec_" + det + ".png", axes, folder);
    createTH2Plot(c1, plot_2D_map["Diff Eta Phi Plots"]["LCs"][det].front(),"Diff_Eta_Phi_CP_LC_" + det + ".png", axes, folder);


    axes = {"#Delta #x", "#Delta #y"};
    createTH2Plot(c1, plot_2D_map["Diff x y Plots"]["Simhits"][det].front(),"Diff_x_y_CP_Sim_" + det + ".png", axes, folder);
    createTH2Plot(c1, plot_2D_map["Diff x y Plots"]["Rechits"][det].front(),"Diff_x_y_CP_Rec_" + det + ".png", axes, folder);
    createTH2Plot(c1, plot_2D_map["Diff x y Plots"]["LCs"][det].front(),"Diff_x_y_CP_LC_" + det + ".png", axes, folder);

    /*

    for(int i = 0; i<47; i++){

      // Layer hits
      std::stringstream ss;
      ss << "Layer" <<std::setw(2) << std::setfill('0') << i <<"_"<<det<<".png";
      t_name = ss.str();

      std::stringstream title;
      title << "Hits in Layer " << std::setfill('0') << i;
      ttitle = title.str();
      hists = {plot_map["Hit Layer Plots"]["Simhits"][det][i],plot_map["Hit Layer Plots"]["Rechits"][det][i], plot_map["Hit Layer Plots"]["LCs"][det][i]}; //i-1 necessary so to map position '0' to 'Layer1'
      createTHSPlot(c1, colors, hists, t_name, axes, {"Simhits", "Rechits", "LCs"}, pos, folder, ttitle);


      // Layer Distance

      ss.str("");
      ss << "_Layer" <<std::setw(2) << std::setfill('0') << i <<"_" << det <<".png";
      t_name = ss.str();

      axes = {"#Delta x", "#Delta y"};
      createTH1Plot(c1,plot_map["Diff x Layer Plots"]["Simhits"][det][i], "Diff_x_CP_Sim"+t_name, {"#Delta x","Occurence"}, folder+"/Layerwise_Distance");
      createTH1Plot(c1,plot_map["Diff y Layer Plots"]["Simhits"][det][i], "Diff_y_CP_Sim"+t_name, {"#Delta y","Occurence"}, folder+"/Layerwise_Distance");
      createTH1Plot(c1,plot_map["Diff Eta Layer Plots"]["Simhits"][det][i], "Diff_eta_CP_Sim"+t_name, {"#Delta #eta","Occurence"}, folder+"/Layerwise_Distance");
      createTH1Plot(c1,plot_map["Diff Phi Layer Plots"]["Simhits"][det][i], "Diff_phi_CP_Sim"+t_name, {"#Delta #phi","Occurence"}, folder+"/Layerwise_Distance");

      createTH1Plot(c1,plot_map["Diff x Layer Plots"]["Rechits"][det][i], "Diff_x_CP_rec"+t_name, {"#Delta x","Occurence"}, folder+"/Layerwise_Distance");
      createTH1Plot(c1,plot_map["Diff y Layer Plots"]["Rechits"][det][i], "Diff_y_CP_rec"+t_name, {"#Delta y","Occurence"}, folder+"/Layerwise_Distance");
      createTH1Plot(c1,plot_map["Diff Eta Layer Plots"]["Rechits"][det][i], "Diff_eta_CP_rec"+t_name, {"#Delta #eta","Occurence"}, folder+"/Layerwise_Distance");
      createTH1Plot(c1,plot_map["Diff Phi Layer Plots"]["Rechits"][det][i], "Diff_phi_CP_rec"+t_name, {"#Delta #phi","Occurence"}, folder+"/Layerwise_Distance");

      createTH1Plot(c1,plot_map["Diff x Layer Plots"]["LCs"][det][i], "Diff_x_CP_lc"+t_name, {"#Delta x","Occurence"}, folder+"/Layerwise_Distance");
      createTH1Plot(c1,plot_map["Diff y Layer Plots"]["LCs"][det][i], "Diff_y_CP_lc"+t_name, {"#Delta y","Occurence"}, folder+"/Layerwise_Distance");
      createTH1Plot(c1,plot_map["Diff Eta Layer Plots"]["LCs"][det][i], "Diff_eta_CP_lc"+t_name, {"#Delta #eta","Occurence"}, folder+"/Layerwise_Distance");
      createTH1Plot(c1,plot_map["Diff Phi Layer Plots"]["LCs"][det][i], "Diff_phi_CP_lc"+t_name, {"#Delta #phi","Occurence"}, folder+"/Layerwise_Distance");
    }

  */
  /*
  }
  */
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

  std::cout << "analyze" << std::endl;

  edm::Handle<std::vector<ticl::Trackster>> ticlTrackstersMerge;
  iEvent.getByToken(ticlTrackstersMergeToken, ticlTrackstersMerge);
  const std::vector<ticl::Trackster>& tracksters = *ticlTrackstersMerge; 


  std::cout<<"Size of TICL Tracksters: " << tracksters.size()<<std::endl;
  std::cout<<"Regressed Energy of Trackster: " << tracksters.front().regressed_energy()<<std::endl;
  for(int i=0;i<7;i++){
    std::cout<<"PID of Trackster: " << tracksters.front().id_probabilities(i)<<std::endl;
  //std::cout<<"Size of TICL Tracksters: " << tracksters[0].trackPtr().size()<<std::endl;
  }



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

  typedef std::map<std::string, std::vector<int>> layerdetmap;
  typedef std::map<std::string, int> detmap;

  typedef std::map<std::string, layerdetmap> layerobjmap;
  typedef std::map<std::string, detmap> objmap;

  std::map<std::string, layerobjmap> counter_map;
  std::map<std::string, objmap> connectivity_counter_map;
  objmap skip_counter_map;
  

  std::vector<std::string> objects = {"Simhits", "Rechits", "LCs"};
  std::vector<std::string> detectors = {"", "Si", "Si 120", "Si 200", "Si 300", "Sc"};

  TString t_name;

  for(const auto& obj : objects){
    for(const auto& det : detectors){
      std::vector<int> counter(47,0);
      counter_map["Hit Counters"][obj][det]=counter;
      connectivity_counter_map["Connectivity"][obj][det]=0;
      connectivity_counter_map["Skip"][obj][det]=0;
      connectivity_counter_map["Miss"][obj][det]=0;
      connectivity_counter_map["Max"][obj][det]=0;
      skip_counter_map[obj][det]=0;
    }
  }

  
  float cp_eta = 0;
  float cp_phi = 0;
  
  // Loop over Caloparticles 


  std::vector<DetId> tmprechits_; tmprechits_.clear();
  for (const auto& it_cp : cps) {
    // do something with track parameters, e.g, plot the charge.
    // int charge = track.charge();
    nCalo++;
    const CaloParticle& cp = ((it_cp)); 
    cp_eta = cp.eta();
    cp_phi = cp.phi();

    const SimClusterRefVector& simclusters = cp.simClusters();
    // std::cout<<"CP number of Simhits:" << simclusters.numberOfSimHits() << std::endl;    

    for (const auto& it_simc : simclusters){
      nSim++;
      const SimCluster& simc = (*(it_simc));
      const auto& sc_haf = simc.hits_and_fractions();


      for (const auto& it_sc_haf : sc_haf){
        nSimhits++;
        DetId detid_ = (it_sc_haf.first);
        std::map<DetId,const HGCRecHit *>::const_iterator itcheck = hitMap.find(detid_);
        unsigned int layer_ = recHitTools_.getLayerWithOffset(detid_);   

        // Investigation of Scintillator-Eta 2.9 discrepancy

        float edist = cp_eta - recHitTools_.getPosition(detid_).eta();
        float phidist = cp_phi - recHitTools_.getPosition(detid_).phi();
        // float drdist = getDr(recHitTools_.getPosition(detid_).eta(),recHitTools_.getPosition(detid_).phi(),cp_eta,cp_phi);

        float cpPositionAtSimHit_x = -1.*recHitTools_.getPosition(detid_).z()*TMath::Tan(2.*(TMath::ATan(TMath::Exp(cp_eta))))*TMath::Cos(cp_phi);
        float cpPositionAtSimHit_y = -1.*recHitTools_.getPosition(detid_).z()*TMath::Tan(2.*(TMath::ATan(TMath::Exp(cp_eta))))*TMath::Sin(cp_phi);

        float xdist = cpPositionAtSimHit_x - recHitTools_.getPosition(detid_).x();
        float ydist = cpPositionAtSimHit_y - recHitTools_.getPosition(detid_).y();

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

        std::string obj = "Simhits";
        for(const auto& det : detectors){
          if(det == ""|| det == detector || det.find(thickness)!=string::npos){
            counter_map["Hit Counters"][obj][det][layer_-1]++;
            plot_map["Diff Eta Plots"][obj][det].front()->Fill(edist);
            plot_map["Diff Phi Plots"][obj][det].front()->Fill(phidist);
            plot_map["Diff x Plots"][obj][det].front()->Fill(xdist);
            plot_map["Diff y Plots"][obj][det].front()->Fill(ydist);
            plot_map["Diff x Layer Plots"][obj][det][layer_-1]->Fill(xdist);
            plot_map["Diff y Layer Plots"][obj][det][layer_-1]->Fill(ydist);
            plot_map["Diff Eta Layer Plots"][obj][det][layer_-1]->Fill(edist);
            plot_map["Diff Phi Layer Plots"][obj][det][layer_-1]->Fill(phidist);

            plot_2D_map["Diff Eta Phi Plots"][obj][det].front()->Fill(edist,phidist); 
            plot_2D_map["Diff x y Plots"][obj][det].front()->Fill(xdist,ydist); 


          }
        }

        std::cout << "Rec" << std::endl;
        if (itcheck != hitMap.end()){
          tmprechits_.push_back(detid_);
          nRechits++;

          // Fill energy histograms

          plot_map["Energy"]["Rechits"][""].front()->Fill((it_sc_haf.second)*itcheck->second->energy());

          // Not necessary? Detid of Simhit and Rechit the same.

          std::string obj = "Rechits";
          for(const auto& det : detectors){
            if(det == "" || det == detector || det.find(thickness)!=string::npos){
              counter_map["Hit Counters"][obj][det][layer_-1]++;

              plot_map["Diff Eta Plots"][obj][det].front()->Fill(edist);
              plot_map["Diff Phi Plots"][obj][det].front()->Fill(phidist);
              plot_map["Diff x Plots"][obj][det].front()->Fill(xdist);
              plot_map["Diff y Plots"][obj][det].front()->Fill(ydist);

              plot_map["Diff x Layer Plots"][obj][det][layer_-1]->Fill(xdist);
              plot_map["Diff y Layer Plots"][obj][det][layer_-1]->Fill(ydist);
              plot_map["Diff Eta Layer Plots"][obj][det][layer_-1]->Fill(edist);
              plot_map["Diff Phi Layer Plots"][obj][det][layer_-1]->Fill(phidist);

              plot_2D_map["Diff Eta Phi Plots"][obj][det].front()->Fill(edist,phidist); 
              plot_2D_map["Diff x y Plots"][obj][det].front()->Fill(xdist,ydist); 

            }
          }
        }
      }
    }
  }
  

  // Loop over LCs  

  std::cout << "LC" << std::endl;
  for (const auto& it_lc : lcs) {

    //if(M_PI < cp_eta && cp_phi <0){
    //  continue;
    //}

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
	      // float drdist = getDr(recHitTools_.getPosition(detid_).eta(),recHitTools_.getPosition(detid_).phi(),cp_eta,cp_phi);

        float cpPositionAtLC_x = -1.*it_lc.z()*TMath::Tan(2.*(TMath::ATan(TMath::Exp(cp_eta))))*TMath::Cos(cp_phi);
        float cpPositionAtLC_y = -1.*it_lc.z()*TMath::Tan(2.*(TMath::ATan(TMath::Exp(cp_eta))))*TMath::Sin(cp_phi);

        float xdist = cpPositionAtLC_x - recHitTools_.getPosition(detid_).x();
        float ydist = cpPositionAtLC_y - recHitTools_.getPosition(detid_).y();

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

        std::string obj = "LCs";
        for(const auto& det : detectors){
          if(det == "" || det == detector || det.find(thickness)!=string::npos){
            std::cout<<det<<"\t"<<detector<<"\t"<<thickness<<std::endl;
            std::cout<<(det=="")<<std::endl;
            std::cout<<(det==detector)<<std::endl;
            std::cout<<(det.find(thickness)!=string::npos)<<std::endl;


            counter_map["Hit Counters"][obj][det][layer_-1]++;

            plot_map["Diff Eta Plots"][obj][det].front()->Fill(edist);
            plot_map["Diff Phi Plots"][obj][det].front()->Fill(phidist);
            plot_map["Diff x Plots"][obj][det].front()->Fill(xdist);
            plot_map["Diff y Plots"][obj][det].front()->Fill(ydist);

            plot_map["Diff x Layer Plots"][obj][det][layer_-1]->Fill(xdist);
            plot_map["Diff y Layer Plots"][obj][det][layer_-1]->Fill(ydist);
            plot_map["Diff Eta Layer Plots"][obj][det][layer_-1]->Fill(edist);
            plot_map["Diff Phi Layer Plots"][obj][det][layer_-1]->Fill(phidist);

            plot_2D_map["Diff Eta Phi Plots"][obj][det].front()->Fill(edist,phidist); 
            plot_2D_map["Diff x y Plots"][obj][det].front()->Fill(xdist,ydist); 
          }
        }
      }
    }
  }
  
  cout<<"nTrack = "<<nTrack<<endl;
  cout<<"nCaloparticles = "<<nCalo<<endl;
  cout<<"nSimparticles = "<<nSim<<endl;
  cout<<"nLCparticles = "<<nLC<<endl;
  cout<<"nSimhits = "<<nSimhits<<endl;
  cout<<"nRechits = "<<nRechits<<endl;
  cout << "Num of LCs = "<< lcs.size() <<endl ;
  cout<<"nLChits = "<<nLChits<<endl;

  //plot_map["Detector Hit Plots"]["Simhits"][""].front()->Fill(nSimhits);
  //plot_map["Detector Hit Plots"]["Rechits"][""].front()->Fill(nRechits);
  //plot_map["Detector Hit Plots"]["LCs"][""].front()->Fill(nLChits);

  plot_map["Validation Plots"]["Simhits"][""].front()->Fill(nSim);
  plot_map["Validation Plots"]["LCs"][""].front()->Fill(nLC);

  for(int i=0; i<47;i++){
    
    // Missing Simhits

    if(counter_map["Hit Counters"]["Simhits"][""][i]==0){
      std::cout << iEvent.id().event() << std::endl;
      missing_simhits_histo->Fill(i,iEvent.id().event());
    }

    for(const auto& obj : objects){
      for(const auto& det : detectors){

        plot_map["Hit Layer Plots"][obj][det][i]->Fill(counter_map["Hit Counters"][obj][det][i]);
        plot_map["Detector Hit Plots"][obj][det].front()->Fill(i,counter_map["Hit Counters"][obj][det][i]);

        // Connectivity Histograms

        //if(hit_layer_counters[j][i]!=0){
        if(counter_map["Hit Counters"][obj][det][i]!=0){

          plot_map["Detector Bool Hit Plots"][obj][det].front()->Fill(i);
          connectivity_counter_map["Connectivity"][obj][det]++;
          connectivity_counter_map["Skip"][obj][det]++;
          skip_counter_map[obj][det]++;
          if(connectivity_counter_map["Miss"][obj][det]!=0){plot_map["Miss Connectivity Plots"][obj][det].front()->Fill(connectivity_counter_map["Miss"][obj][det]);}
          connectivity_counter_map["Miss"][obj][det]=0;
          if(i==46){
            plot_map["Connectivity Plots"][obj][det].front()->Fill(connectivity_counter_map["Connectivity"][obj][det]);
            plot_map["Skip Connectivity Plots"][obj][det].front()->Fill(connectivity_counter_map["Skip"][obj][det]);
            if(connectivity_counter_map["Max"][obj][det]<connectivity_counter_map["Connectivity"][obj][det]){
              connectivity_counter_map["Max"][obj][det] = connectivity_counter_map["Connectivity"][obj][det];
            }
            plot_map["Max Connectivity Plots"][obj][det].front()->Fill(connectivity_counter_map["Max"][obj][det]);


            for(int k=i+1-connectivity_counter_map["Connectivity"][obj][det]; k<i+1; k++){
              plot_map["W Connectivity Plots"][obj][det].front()->Fill(k, connectivity_counter_map["Connectivity"][obj][det]);
            }
          }
        }

        else{
          connectivity_counter_map["Miss"][obj][det]++;
          if(skip_counter_map[obj][det]!=0){
            skip_counter_map[obj][det]--;
            connectivity_counter_map["Skip"][obj][det]++;
          }
          else{
            if(connectivity_counter_map["Skip"][obj][det]-skip_>0){plot_map["Skip Connectivity Plots"][obj][det].front()->Fill(connectivity_counter_map["Connectivity"][obj][det]);}

            connectivity_counter_map["Skip"][obj][det]=0;
          }


          if(connectivity_counter_map["Connectivity"][obj][det]!=0){
            plot_map["Connectivity Plots"][obj][det].front()->Fill(connectivity_counter_map["Connectivity"][obj][det]);
          }

          for(int k=i-connectivity_counter_map["Connectivity"][obj][det]; k<i; k++){
            plot_map["W Connectivity Plots"][obj][det].front()->Fill(k,connectivity_counter_map["Connectivity"][obj][det]);
          }
          if(connectivity_counter_map["Max"][obj][det]<connectivity_counter_map["Connectivity"][obj][det]){
            connectivity_counter_map["Max"][obj][det]=connectivity_counter_map["Connectivity"][obj][det];
          }
          connectivity_counter_map["Connectivity"][obj][det]=0;
          if(i==46){
            plot_map["Max Connectivity Plots"][obj][det].front()->Fill(connectivity_counter_map["Max"][obj][det]);
            plot_map["Miss Connectivity Plots"][obj][det].front()->Fill(connectivity_counter_map["Miss"][obj][det]);
          }
        }
      }
    }
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
