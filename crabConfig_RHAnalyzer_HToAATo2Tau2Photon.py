from CRABClient.UserUtilities import config#, getUsernameFromSiteDB
config = config()
# See parameter defintions here: https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile#CRAB_configuration_parameters
# To test crab config:
# crab submit -c crabConfig_data.py --dryrun
# To submit to crab:
# crab submit -c crabConfig_data.py
# To check job status:
# crab status -d <config.General.workArea>/<config.General.requestName>
# To resubmit jobs:
# crab resubmit -d <config.General.workArea>/<config.General.requestName>
Mass_tag = 'm3p6To8'
# Local job directory will be created in:
# <config.General.workArea>/<config.General.requestName>
config.General.workArea        = 'large_scale_runs'
config.General.requestName     = 'HToAATo2TauHad2Photon_RHAnalyzer_%s'%Mass_tag
config.General.transferOutputs = True
config.General.transferLogs    = True

# CMS cfg file goes here:
config.JobType.pluginName  = 'Analysis' # mass > 8 use this
config.JobType.psetName    = 'RecHitAnalyzer/python/ConfFile_cfg.py' # cms cfg file for generating events
# config.JobType.maxMemoryMB = 5000 #5000


config.JobType.maxMemoryMB = 4000
config.JobType.numCores = 8

config.Data.inputDBS = 'phys03'
config.JobType.allowUndistributedCMSSW = True
# Define input and units per job here:
dataset  = '/GEN_SIM_HToAATo2Tau2Photon_m3p6To8_m3p6To8/lpcml-m3p6To8_MiniAODSIM-d639958dfa8e672fc82464be0b06406f/USER'


config.Data.inputDataset   = dataset
config.Data.splitting      = 'FileBased'
config.Data.unitsPerJob    = 1  # units: as defined by config.Data.splitting
config.Data.totalUnits     = -1 # -1: all inputs. total jobs submitted = totalUnits / unitsPerJob. cap of 10k jobs per submission
config.Data.publication    = False



# Output files will be stored in config.Site.storageSite at directory:
# <config.Data.outLFNDirBase>/<config.Data.outputPrimaryDataset>/<config.Data.outputDatasetTag>/
config.Data.outLFNDirBase = '/store/user/ereinhar/MCGeneration_run3'

config.Site.storageSite = 'T3_US_FNALLPC'
config.Data.outputDatasetTag = config.General.requestName

config.Data.ignoreLocality = True
config.Site.whitelist = [
    'T2_AT_Vienna', 'T2_BE_IIHE', 'T2_BE_UCL', 'T2_BR_SPRACE', 'T2_BR_UERJ',
    'T2_CH_CERN', 'T2_CN_Beijing', 'T2_DE_DESY', 'T2_DE_RWTH',
    'T2_EE_Estonia', 'T2_ES_CIEMAT', 'T2_ES_IFCA', 'T2_FI_HIP', 
    'T2_FR_IPHC', 'T2_GR_Ioannina', 'T2_HU_Budapest', 'T2_IN_TIFR',
    'T2_IT_Bari', 'T2_IT_Legnaro', 'T2_IT_Pisa', 'T2_IT_Rome',
    'T2_KR_KISTI', 'T2_PK_NCP', 'T2_PL_Cyfronet', 
    'T2_PT_NCG_Lisbon', 'T2_RU_IHEP', 
    'T2_TR_METU', 'T2_TW_NCHC', 'T2_UA_KIPT',
    'T2_UK_London_Brunel', 'T2_UK_London_IC', 'T2_UK_SGrid_Bristol',
    'T2_UK_SGrid_RALPP', 'T2_US_Caltech', 'T2_US_Florida',
    'T2_US_MIT', 'T2_US_Nebraska', 'T2_US_Purdue', 'T2_US_UCSD',
    'T2_US_Vanderbilt', 'T2_US_Wisconsin'
]