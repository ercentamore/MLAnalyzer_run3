import os

#cfg='RecHitAnalyzer/python/ConfFile_data_cfg.py'
cfg='RecHitAnalyzer/python/ConfFile_cfg.py'
inputFiles_='file:/uscms/home/ecentamo/nobackup/CMSSW_13_0_13/src/nobackup/CMSSW_13_0_13/src/MCProduction/E2E-HToAATo2Tau2Photon/NanoAOD_HToAATo2Gluon2Photon.root'#pixel checks
# inputFiles_='file:step3_AODSIM_M14_1.root'

maxEvents_=-1
skipEvents_=0#
# outputFile_='MLAnalyzer_HToAATo2Tau2Photon_with_triggers.root'
outputFile_='MLAnalyzer_HToAATo2Gluon2Photon.root'

# cmd="cmsTraceExceptions cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d outputFile=%s"%(cfg,inputFiles_,maxEvents_,skipEvents_,outputFile_)
cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d outputFile=%s"%(cfg,inputFiles_,maxEvents_,skipEvents_,outputFile_)
print(f"{cmd}")
os.system(cmd)
