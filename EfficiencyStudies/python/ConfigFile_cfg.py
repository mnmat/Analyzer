import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
process = cms.Process('PROD',Phase2C17I13M9)
process.load('Configuration.Geometry.GeometryExtended2026D86Reco_cff')


process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
            'file:/afs/cern.ch/user/m/mmatthew/CMSSW_12_3_0_pre5/src/hgcal_private/production/None/singlemuon_flatEGun_hgcalCenter/step3/step3_singlemuon_e100GeV_nopu.root'
                )
                            )

process.demo = cms.EDAnalyzer('EfficiencyStudies',
   tracks = cms.untracked.InputTag('generalTracks'),
   caloParticles = cms.InputTag("mix", "MergedCaloTruth"),
   hgcalRecHitsEE = cms.InputTag("HGCalRecHit", "HGCEERecHits"),
   hgcalRecHitsFH = cms.InputTag("HGCalRecHit", "HGCHEFRecHits"),
   hgcalRecHitsBH = cms.InputTag("HGCalRecHit", "HGCHEBRecHits"),
   hgcalLayerClusters = cms.InputTag("hgcalLayerClusters", "", "RECO"),
   trackPtMin = cms.double(0.3)
                              )

process.p = cms.Path(process.demo)
