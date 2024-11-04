import os

#cfg='RecHitAnalyzer/python/ConfFile_data_cfg.py'
cfg='RecHitAnalyzer/python/ConfFile_cfg.py'
inputFiles_='file:/uscms/home/bbbam/nobackup/analysis_run3/MCGeneration/CMSSW_13_0_17/src/MCProduction_run3/E2E-ATo2Tau/AOD_ATo2Tau_extra_collection.root'
# inputFiles_='root://cmseos.fnal.gov//store/group/lpcml/bbbam/MCGeneration/signal_withTrigger/GEN_SIM_HToAATo4Tau_M_8_pythia8_2018UL/raw_to_AODSIM_HToAATo4Tau_M_8/240606_182926/0000/AODSIM_HToAATo4Tau_46_6.root'#pixel checks

#maxEvents_=10
#maxEvents_=20
maxEvents_=-1
skipEvents_=0#
#outputFile_='MLAnal_PhaseI_TTbar_13TeVu_trackRefitter.root'
#outputFile_='GJet.root'
#outputFile_='ttbar_secVertex.root'
#outputFile_='DYToTauTau_subJet.root'
#outputFile_='WJets_secVertex.root'
#outputFile_='dyToEE.root'
#outputFile_='acd_EmEnriched.root'
outputFile_='Ato2Tau_massreg_sample.root'

# cmd="cmsTraceExceptions cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d outputFile=%s"%(cfg,inputFiles_,maxEvents_,skipEvents_,outputFile_)
cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d outputFile=%s"%(cfg,inputFiles_,maxEvents_,skipEvents_,outputFile_)
print(f"{cmd}")
os.system(cmd)
