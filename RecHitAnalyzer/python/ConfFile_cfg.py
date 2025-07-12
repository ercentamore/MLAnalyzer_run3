
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('analysis')
options.register('skipEvents',
    default=0,
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.int,
    info = "skipEvents")
# TODO: put this option in cmsRun scripts
options.register('processMode',
    default='JetLevel',
    #default='EventLevel',
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.string,
    info = "process mode: JetLevel or EventLevel")
options.parseArguments()


from Configuration.Eras.Era_Run3_cff import Run3
process = cms.Process("FEVTAnalyzer", Run3)
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
process.load("RecoLocalTracker.SiPixelRecHits.SiPixelRecHits_cfi")
process.load("RecoLocalTracker.SiStripRecHitConverter.SiStripRecHitConverter_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('130X_mcRun3_2023_realistic_postBPix_v5')
# process.GlobalTag.globaltag = cms.string('130X_dataRun3_HLT_v1')
process.es_prefer_GlobalTag = cms.ESPrefer('PoolDBESSource','GlobalTag')

process.load("RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi")

process.load("JetMETCorrections.Modules.JetResolutionESProducer_cfi")
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")

process.TrackRefitter.TTRHBuilder = 'WithAngleAndTemplate'
process.MeasurementTrackerEvent.inactivePixelDetectorLabels = cms.VInputTag() # to prebent worning:  fail to get the list of inactive pixel modules, because of 4.2/4.4 event content change. 


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
    )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      options.inputFiles
      )
    , skipEvents = cms.untracked.uint32(options.skipEvents)
    )
print (" >> Loaded",len(options.inputFiles),"input files from list.")


process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
    # ,SkipEvent = cms.untracked.vstring('ProductNotFound')
)

process.load("MLAnalyzer_run3.RecHitAnalyzer.RHAnalyzer_cfi")
process.fevt.mode = cms.string(options.processMode)
print (" >> Processing as:",(process.fevt.mode))

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outputFile)
    )

############################
# Event Analysis
############################

process.hltFilter = cms.EDFilter("HLTHighLevel",
                                eventSetupPathsKey = cms.string(''),
                                TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
                                # HLTPaths = cms.vstring('*'),
                                HLTPaths = cms.vstring('HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_v','HLT_DoublePhoton33_CaloIdL_v'),
                                andOr = cms.bool(True),
                                throw = cms.bool(False)
                                )

process.p = cms.Path(
process.siStripMatchedRecHits*process.siPixelRecHits*process.MeasurementTrackerEvent*process.TrackRefitter*
#process.hltFilter*
process.fevt
)
