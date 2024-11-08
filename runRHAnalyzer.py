import os

#cfg='RecHitAnalyzer/python/ConfFile_data_cfg.py'
cfg='RecHitAnalyzer/python/ConfFile_cfg.py'
# inputFiles_='file:/uscms/home/bbbam/nobackup/analysis_run3/MCGeneration/CMSSW_13_0_17/src/MCProduction_run3/E2E-ATo2Tau/AOD_ATo2Tau_extra_collection.root'
# inputFiles_='root://cmseos.fnal.gov//store/group/lpcml/bbbam/MCGeneration_run3/GEN_SIM_ATo2Tau_m1p2To3p6_pt30To300_v4/AOD_ATo4Tau_Hadronic_m1p2To3p6/241103_223217/0001/AOD_ATo2Tau_extra_collection_1916.root'#pixel checks
inputFiles_='root://cmseos.fnal.gov//store/group/lpcml/bbbam/MCGeneration_run3/GEN_SIM_ATo2Tau_m3p6To18_pt30To300_v2/AOD_ATo4Tau_Hadronic_m3p6To18/241103_222918/0000/AOD_ATo2Tau_extra_collection_993.root'#pixel checks

#maxEvents_=10
#maxEvents_=20
maxEvents_=2
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
