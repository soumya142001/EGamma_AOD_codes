import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
                               #'file:/eos/user/s/sosarkar/1EE5922B-BCBD-2844-82A1-2A6BB9B4ACD3.root'
                               #'file:/eos/user/s/sosarkar/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_AOD.root'
                              'file:/eos/user/s/sosarkar/BuToKJPsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV_pythia8_Run3_Summer_2021_120x_AODSIM.root',
                              'file:/eos/user/s/sosarkar/BuToKJPsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV_pythia8_Run3_Summer_2021_120x_AODSIM_2.root',
                              'file:/eos/user/s/sosarkar/BuToKJPsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV_pythia8_Run3_Summer_2021_120x_AODSIM_3.root',
                              'file:/eos/user/s/sosarkar/BuToKJPsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV_pythia8_Run3_Summer_2021_120x_AODSIM_4.root',
                              'file:/eos/user/s/sosarkar/BuToKJPsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV_pythia8_Run3_Summer_2021_120x_AODSIM_5.root',
                              'file:/eos/user/s/sosarkar/BuToKJPsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV_pythia8_Run3_Summer_2021_120x_AODSIM_6.root',
                              'file:/eos/user/s/sosarkar/BuToKJPsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV_pythia8_Run3_Summer_2021_120x_AODSIM_7.root'
                               # 'file:/eos/user/s/sosarkar/DYJetsToLL_M-50_TuneCP5_RunIISummer20UL18RECO_amcantnloFXFX-pythia8_AOD.root'
                )
                            )

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('EGamma_JPsi_analysis_v2.root')
                                   )

process.demo = cms.EDAnalyzer('Analyzer_egamma',
   
   gsf_electron_ged = cms.InputTag('gedGsfElectrons'),
   gsf_electron_unonly = cms.InputTag('uncleanedOnlyGsfElectrons'),
   conversions = cms.InputTag('allConversions'),
   DeDx_token = cms.InputTag('dedxHitInfo'),
   track_general = cms.InputTag('generalTracks'),
   track_displaced = cms.InputTag('displacedTracks'),
   GsfTrack_token = cms.InputTag('electronGsfTracks'),
   sc_token = cms.InputTag('particleFlowEGamma'),
   calo_token = cms.InputTag('particleFlowEGamma','EBEEClusters'),
   calo_eb_token = cms.InputTag('particleFlowSuperClusterECAL','particleFlowBasicClusterECALBarrel'),
   calo_ee_token = cms.InputTag('particleFlowSuperClusterECAL','particleFlowBasicClusterECALEndcap'),
   ecalrechit = cms.InputTag('reducedEcalRecHitsEB'),
   photon_token = cms.InputTag('gedPhotons'),
   gen_token = cms.InputTag('genParticles')                           
   #tracks    = cms.untracked.InputTag('generalTracks'),
   #trackPtMin = cms.double(0.3),
   #trackEtaMin = cms.double(-2.4),
   #trackEtaMax = cms.double(2.4)
                              )

process.p = cms.Path(process.demo)
