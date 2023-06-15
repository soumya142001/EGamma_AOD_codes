from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config
config = Configuration()

config.section_("General")
config.General.requestName      = 'DYJetsToLL_M50_UL2018'
config.General.transferLogs     = True

config.section_("JobType")
config.JobType.allowUndistributedCMSSW = True
config.JobType.pluginName       = 'Analysis'
config.JobType.psetName         = '/afs/cern.ch/user/s/sosarkar/MYDEMOANALYZER/CMSSW_12_0_0/src/Demo_e_gamma/Analyzer_egamma/python/ConfFile_cfg.py'
config.JobType.sendPythonFolder = True

config.section_("Data")
config.Data.inputDataset        = '/DYJetsToLL_M-500to700_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v3/AODSIM'
config.Data.inputDBS            = 'global'
config.Data.allowNonValidInputDataset = True
config.Data.splitting           = 'FileBased'
config.Data.totalUnits          = 51
config.Data.unitsPerJob         = 5
config.Data.outLFNDirBase       = '/store/user/sosarkar/DYJetsToLL_M50'
config.Data.publication = False
config.Data.outputDatasetTag    = 'DYJetsToLL_M50'

config.section_("Site")
config.Site.storageSite         = "T3_CH_CERNBOX"


