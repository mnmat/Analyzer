import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
process = cms.Process('PROD',Phase2C17I13M9)
process.load('Configuration.Geometry.GeometryExtended2026D86Reco_cff')


process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Specify Eta

eta = "Eta_17_27" # Eta_16, Eta_17, Eta_29, Eta_17_27

"""

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
            'file:/afs/cern.ch/user/m/mmatthew/CMSSW_12_3_0_pre5/src/hgcal_private/production/'+eta+'/singleel_flatEGun_hgcalCenter/step3/step3_singleel_e100GeV_nopu.root'
                )
                            )
outfile_  = 'file:file:/afs/cern.ch/user/m/mmatthew/CMSSW_12_3_0_pre5/src/hgcal_private/production/'+eta+'/singleel_flatEGun_hgcalCenter/step3/ttree_singleel_e100GeV_nopu.root'

"""

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
            'file:/eos/home-m/mmatthew/Data/Analyzer/EfficiencyStudies/Test/' + eta + '/singlemuon_flatEGun_hgcalCenter/step3/step3_singlemuon_e100GeV_nopu.root'
                )
                            )
outfile_  = 'file:file:/eos/home-m/mmatthew/Data/Analyzer/EfficiencyStudies/Test/'+eta+'/singlemuon_flatEGun_hgcalCenter/step3/ttree_singlemuon_e100GeV_nopu.root'


"""
process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
            'file:/afs/cern.ch/user/m/mmatthew/CMSSW_12_3_0_pre5/src/hgcal_private/production/singleel/singleel_flatEGun_hgcalCenter/step3/step3_singleel_e100GeV_nopu.root'
                )
                            )
outfile_  = 'file:file:/afs/cern.ch/user/m/mmatthew/CMSSW_12_3_0_pre5/src/hgcal_private/production/singleel/singleel_flatEGun_hgcalCenter/step3/ttree_singleel_e100GeV_nopu.root'

"""
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(outfile_),
                                   closeFileFast = cms.untracked.bool(True)
                               )


process.demo = cms.EDAnalyzer('EfficiencyStudies',
   tracks = cms.untracked.InputTag('generalTracks'),
   caloParticles = cms.InputTag("mix", "MergedCaloTruth"),
   Tracksters = cms.InputTag("ticlTrackstersMerge","","RECO"),
   hgcalRecHitsEE = cms.InputTag("HGCalRecHit", "HGCEERecHits"),
   hgcalRecHitsFH = cms.InputTag("HGCalRecHit", "HGCHEFRecHits"),
   hgcalRecHitsBH = cms.InputTag("HGCalRecHit", "HGCHEBRecHits"),
   hgcalLayerClusters = cms.InputTag("hgcalLayerClusters", "", "RECO"),
   skip = cms.int32(1),
   eta = cms.string(eta),  
   #trackPtMin = cms.double(0.3)
                              )

process.p = cms.Path(process.demo)
