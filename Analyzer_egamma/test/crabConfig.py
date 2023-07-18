from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config
config = Configuration()

config.section_("General")
config.General.requestName      = 'JPsiToee_04July_v2' 
config.General.transferLogs     = True

config.section_("JobType")
config.JobType.allowUndistributedCMSSW = True
config.JobType.pluginName       = 'Analysis'
config.JobType.psetName         = '/afs/cern.ch/user/s/sosarkar/MYDEMOANALYZER/CMSSW_12_0_0/src/Demo_e_gamma/Analyzer_egamma/python/ConfFile_cfg.py'
config.JobType.sendPythonFolder = True
config.JobType.numCores = 1

config.section_("Data")
config.Data.inputDataset        ='/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/Run3Summer21DR-FlatPU0to70_120X_mcRun3_2021_realistic_v6-v1/AODSIM'
# '/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/RunIIAutumn18RECOBParking-Custom_RK_BParking_for_RK_102X_upgrade2018_realistic_v15-v2/AODSIM'
 # '/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM'
#'/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13p6TeV-pythia8-evtgen/Run3Winter22DR-FlatPU0to70_122X_mcRun3_2021_realistic_v9-v2/AODSIM'
config.Data.inputDBS            = 'global'
config.Data.allowNonValidInputDataset = True
config.Data.splitting           = 'FileBased'
config.Data.totalUnits          = 1
config.Data.unitsPerJob         = 1
config.Data.outLFNDirBase       = '/store/user/sosarkar/JPsiToee_04July_v2'
config.Data.publication = False
config.Data.outputDatasetTag    = 'JPsiToee_04July_v2'

config.section_("Site")
config.Site.storageSite         = "T3_CH_CERNBOX"
config.Site.blacklist           = ['T2_IN_TIFR']
#config.Site.whitelist           = ['T2_TR_METU']
