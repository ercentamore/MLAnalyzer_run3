import os

#cfg='RecHitAnalyzer/python/ConfFile_data_cfg.py'
cfg='RecHitAnalyzer/python/ConfFile_cfg.py'
inputFiles_='file:/uscms/home/ereinhar/nobackup/CMSSW_13_0_13/src/GEN_SIM_HToAATo2Tau2Photon.root'#pixel checks
# inputFiles_='file:step3_AODSIM_M14_1.root'

maxEvents_=-1
skipEvents_=0#
outputFile_='HToAATo2Tau2Photon_MLAnalyzer.root'

# cmd="cmsTraceExceptions cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d outputFile=%s"%(cfg,inputFiles_,maxEvents_,skipEvents_,outputFile_)
cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d outputFile=%s"%(cfg,inputFiles_,maxEvents_,skipEvents_,outputFile_)
print(f"{cmd}")
os.system(cmd)
