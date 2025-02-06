import os

#cfg='RecHitAnalyzer/python/ConfFile_data_cfg.py'
cfg='RecHitAnalyzer/python/ConfFile_cfg.py'
# inputFiles_='file:/uscms/home/bbbam/nobackup/analysis_run3/MCGeneration/CMSSW_13_0_17/src/MCProduction_run3/E2E-ATo2Tau/AOD_AToTau_m1p8To3p6_pt30To300.root'
# inputFiles_='root://cmseos.fnal.gov//store/group/lpcml/bbbam/MCGeneration_run3/GEN_SIM_ATo2Tau_m1p2To3p6_pt30To300_v4/AOD_ATo4Tau_Hadronic_m1p2To3p6/241103_223217/0001/AOD_ATo2Tau_extra_collection_1916.root'#pixel checks
inputFiles_='root://cmseos.fnal.gov//store/group/lpcml/bbbam/MCGeneration_run3/GEN_SIM_Tau_m1p8To3p6_pt30To300/AOD_Tau_Hadronic_decay_m1p8To3p6_v2/250131_023511/0000/AOD_AToTau_m1p8To3p6_pt30To300_88.root'#pixel checks
# inputFiles_='file:step3_AODSIM_M14_1.root'

#maxEvents_=10
maxEvents_=20
# maxEvents_=-1
skipEvents_=0#
#outputFile_='MLAnal_PhaseI_TTbar_13TeVu_trackRefitter.root'
#outputFile_='GJet.root'
#outputFile_='ttbar_secVertex.root'
#outputFile_='DYToTauTau_subJet.root'
#outputFile_='WJets_secVertex.root'
#outputFile_='dyToEE.root'
# outputFile_='QCD_pt15to7000_Run3Summer23GS.root'
outputFile_='AToTau_massreg_sample_m1p8To3p6_pt30T0300.root'
# outputFile_='HToAAto4Tau_signal_sample_m3p7.root'

# cmd="cmsTraceExceptions cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d outputFile=%s"%(cfg,inputFiles_,maxEvents_,skipEvents_,outputFile_)
cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d outputFile=%s"%(cfg,inputFiles_,maxEvents_,skipEvents_,outputFile_)
print(f"{cmd}")
os.system(cmd)
