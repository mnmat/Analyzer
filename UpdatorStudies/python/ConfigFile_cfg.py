import FWCore.ParameterSet.Config as cms
import argparse

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
process = cms.Process('PROD',Phase2C17I13M9)
process.load('Configuration.Geometry.GeometryExtended2026D99Reco_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Specify Eta

import sys

eta = sys.argv[2]
energy = sys.argv[3]
#eta = "16"
#energy = "10"
propagator = "Analytical" #Analytical, RungeKutta
mb = "mb_ngun"
nevents = "500"
cap = "zpos"
#fname = '/eos/home-m/mmatthew/Data/KF/MissingSimhits/CMSSW_13_1_0_pre1/' + cap+'/n'+nevents+'/Eta_'+eta+'/singlemuon_flatEGun_hgcalCenter/step3/'
fname = '/eos/user/m/mmatthew/Data/KF/CMSSW_13_1_0_pre1/Analyzer/UpdatorStudies/' + cap+'/n'+nevents+'/Eta_'+eta+'/singlemuon_flatEGun_hgcalCenter/step3/'


process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring('file:'+fname + 'step3_singlemuon_e'+energy+'GeV_eta'+eta+'_'+cap+'_events'+nevents+'_nopu.root'))

outfile_  = 'file:'+fname + 'ttree_singlemuon_e'+energy+'GeV_eta'+eta+'_'+cap+'_events'+nevents+'_nopu.root'

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(outfile_),
                                   closeFileFast = cms.untracked.bool(True)
                               )

process.demo = cms.EDAnalyzer('UpdatorStudies',
   caloParticles = cms.InputTag("mix", "MergedCaloTruth"),
   Tracksters = cms.InputTag("ticlTrackstersMerge","","RECO"),
   hgcalRecHitsEE = cms.InputTag("HGCalRecHit", "HGCEERecHits"),
   hgcalRecHitsFH = cms.InputTag("HGCalRecHit", "HGCHEFRecHits"),
   hgcalRecHitsBH = cms.InputTag("HGCalRecHit", "HGCHEBRecHits"),
   KFHits = cms.InputTag("ticlTrackstersKF","KFHits","RECO"),
   PropHits = cms.InputTag("ticlTrackstersKF","PropHits","RECO"),
   hgcalLayerClusters = cms.InputTag("hgcalLayerClusters", "", "RECO"),
   eta = cms.string(eta),
   energy = cms.string(energy), 
   outdir = cms.string(fname), 
   #trackPtMin = cms.double(0.3)
                              )
process.p = cms.Path(process.demo)

