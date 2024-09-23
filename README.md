### Instruction to set up E2E-Analayzer for Run 3
```
export SCRAM_ARCH=el8_amd64_gcc11

source /cvmfs/cms.cern.ch/cmsset_default.sh

cmsrel CMSSW_13_0_13

cd CMSSW_13_0_13/src

eval `scram runtime -sh`

cmsenv
```
# Cloning this repo is neccessary
git clone https://github.com/svfit/ClassicSVfit TauAnalysis/ClassicSVfit -b fastMTT_21_06_2018

git clone https://github.com/svfit/SVfitTF TauAnalysis/SVfitTF

# Clone Actual Analyzer
git clone git@github.com:bhbam/MLAnalyzer_run3.git

scram b -j 16

cd MLAnalyzer_run3

python3 runRHAnalyzer.py # this equivalent to cmsRun command (pythonic version only)
