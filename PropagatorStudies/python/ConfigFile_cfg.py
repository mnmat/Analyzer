import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
process = cms.Process('PROD',Phase2C17I13M9)
process.load('Configuration.Geometry.GeometryExtended2026D86Reco_cff')

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Specify Eta

eta = "29" # Eta_16, Eta_17, Eta_29, Eta_17_27
energy = "100"
propagator = "Analytical" #Analytical, RungeKutta
mb = "10"
fname = 'file:/eos/home-m/mmatthew/Data/Analyzer/PropagatorStudies/HGCTracker/Eta'+eta+'/singlemuon_flatEGun_hgcalCenter/step3/'
#fname = 'file:/eos/home-m/mmatthew/Data/Analyzer/PropagatorStudies/'+propagator+'mb_AG/singlemuon_flatEGun_hgcalCenter/step3/'


process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(fname + 'step3_singlemuon_e'+energy+'GeV_eta'+eta+'_nopu.root'))



#outfile_  = 'file:'+fname + 'ttree_singlemuon_e'+energy+'GeV_eta'+eta+'_nopu.root'
outfile_ = 'file:/eos/home-m/mmatthew/Data/deleteme.root'

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(outfile_),
                                   closeFileFast = cms.untracked.bool(True)
                               )



process.demo = cms.EDAnalyzer('PropagatorStudies',
   tracks = cms.untracked.InputTag('generalTracks'),
   caloParticles = cms.InputTag("mix", "MergedCaloTruth"),
   Tracksters = cms.InputTag("ticlTrackstersMerge","","RECO"),
   hgcalRecHitsEE = cms.InputTag("HGCalRecHit", "HGCEERecHits"),
   hgcalRecHitsFH = cms.InputTag("HGCalRecHit", "HGCHEFRecHits"),
   hgcalRecHitsBH = cms.InputTag("HGCalRecHit", "HGCHEBRecHits"),
   #propagatorEM = cms.InputTag("ticlTrackstersEM","TEST","RECO"),
   #propagatorHAD = cms.InputTag("ticlTrackstersHAD","TEST","RECO"),
   propagatorKF = cms.InputTag("ticlTrackstersKF","TEST","RECO"),
   #propagatorTrk = cms.InputTag("ticlTrackstersTrk","TEST","RECO"),
   #propagatorTrkEM = cms.InputTag("ticlTrackstersTrkEM","TEST","RECO"),
   hgcalLayerClusters = cms.InputTag("hgcalLayerClusters", "", "RECO"),
   eta = cms.string(eta),  
   #trackPtMin = cms.double(0.3)
                              )
process.p = cms.Path(process.demo)
