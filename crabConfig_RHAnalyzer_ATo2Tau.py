from CRABClient.UserUtilities import config#, getUsernameFromSiteDB
config = config()
# See parameter defintions here: https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile#CRAB_configuration_parameters
# To submit to crab:
# crab submit -c crabConfig_data.py
# To check job status:
# crab status -d <config.General.workArea>/<config.General.requestName>
# To resubmit jobs:
# crab resubmit -d <config.General.workArea>/<config.General.requestName>
Mass_tag = 'm1p2To3p6' #
# Local job directory will be created in:
# <config.General.workArea>/<config.General.requestName>
config.General.workArea        = 'crab_RHAnalyzer'
config.General.requestName     = 'RHAnalyzer_ATo4Tau_Hadronic_%s'%Mass_tag
config.General.transferOutputs = True
config.General.transferLogs    = True

# CMS cfg file goes here:
config.JobType.pluginName  = 'Analysis' # mass > 8 use this
config.JobType.psetName    = 'RecHitAnalyzer/python/ConfFile_cfg.py' # cms cfg file for generating events
# config.JobType.maxMemoryMB = 5000 #5000


config.Data.inputDBS = 'phys03'
config.JobType.allowUndistributedCMSSW = True
# Define input and units per job here:
dataset  = {
'm3p6To18':'/GEN_SIM_ATo2Tau_m3p6To18_pt30To300_v2/lpcml-AOD_ATo4Tau_Hadronic_m3p6To18-953b1873547799e513f8a43f2c57e3b2/USER'
,'m1p2To3p6':'/GEN_SIM_ATo2Tau_m1p2To3p6_pt30To300_v4/lpcml-AOD_ATo4Tau_Hadronic_m1p2To3p6-953b1873547799e513f8a43f2c57e3b2/USER'
}.get(Mass_tag, None)


config.Data.inputDataset = dataset
config.Data.splitting      = 'FileBased'
config.Data.unitsPerJob    = 1  # units: as defined by config.Data.splitting
config.Data.totalUnits     = -1 # -1: all inputs. total jobs submitted = totalUnits / unitsPerJob. cap of 10k jobs per submission
config.Data.publication    = False



# Output files will be stored in config.Site.storageSite at directory:
# <config.Data.outLFNDirBase>/<config.Data.outputPrimaryDataset>/<config.Data.outputDatasetTag>/
config.Data.outLFNDirBase = '/store/group/lpcml/bbbam/Ntuples_run3'
# config.Data.outLFNDirBase = '/store/user/bhbam/MCGeneration'

config.Site.storageSite = 'T3_US_FNALLPC'
config.Data.outputDatasetTag = config.General.requestName
