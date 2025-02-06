from CRABClient.UserUtilities import config#, getUsernameFromSiteDB
config = config()
# See parameter defintions here: https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile#CRAB_configuration_parameters
# To submit to crab:
# crab submit -c crabConfig_data.py
# To check job status:
# crab status -d <config.General.workArea>/<config.General.requestName>
# To resubmit jobs:
# crab resubmit -d <config.General.workArea>/<config.General.requestName>
Mass_tag = '14'#'3p7','14' #
# Local job directory will be created in:
# <config.General.workArea>/<config.General.requestName>
config.General.workArea        = 'crab_RHAnalyzer'
config.General.requestName     = 'RHAnalyzer_HToAATo4Tau_Hadronic_M%s'%Mass_tag
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
'3p7':'/HToAATo4Tau_hadronic_tauDecay_M3p7_Run3_2023/phys_diffraction-3p7_AODSIM_multiThreads_8Gb_new-953b1873547799e513f8a43f2c57e3b2/USER'
,'14':'/HToAATo4Tau_hadronic_tauDecay_M14_Run3_2023/phys_diffraction-14_AODSIM_multiThreads_8Gb_new-953b1873547799e513f8a43f2c57e3b2/USER'
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
