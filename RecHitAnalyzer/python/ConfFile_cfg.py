
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

process = cms.Process("FEVTAnalyzer")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
process.load("RecoLocalTracker.SiPixelRecHits.SiPixelRecHits_cfi")
process.load("RecoLocalTracker.SiStripRecHitConverter.SiStripRecHitConverter_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

# process.load("Geometry.ForwardGeometry.CastorGeometry_cfi")
# process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
# process.load("JetMETCorrections.Modules.JetResolutionESProducer_cfi")
# process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
# process.load("Geometry.CaloEventSetup.CaloGeometry_cff")
# process.load("Geometry.ForwardGeometry.ForwardGeometry_cfi")
# process.es_prefer_ZdcGeometryFromDBEP = cms.ESPrefer("ZdcGeometryFromDBEP")
# process.es_prefer_HcalHardcodeGeometryEP = cms.ESPrefer("HcalHardcodeGeometryEP")
# process.es_prefer_CastorGeometryFromDBEP = cms.ESPrefer("CastorGeometryFromDBEP")
# process.es_prefer_CaloTowerGeometryFromDBEP = cms.ESPrefer("CaloTowerGeometryFromDBEP")
# process.es_prefer_EcalBarrelGeometryFromDBEP = cms.ESPrefer("EcalBarrelGeometryFromDBEP")
# process.es_prefer_EcalEndcapGeometryFromDBEP = cms.ESPrefer("EcalEndcapGeometryFromDBEP")
# process.es_prefer_EcalPreshowerGeometryFromDBEP = cms.ESPrefer("EcalPreshowerGeometryFromDBEP")

process.GlobalTag.globaltag = cms.string('130X_mcRun3_2023_realistic_postBPix_v5')
# process.GlobalTag.globaltag = cms.string('130X_dataRun3_HLT_v1')
process.es_prefer_GlobalTag = cms.ESPrefer('PoolDBESSource','GlobalTag')

process.TrackRefitter.TTRHBuilder = 'WithAngleAndTemplate'

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

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

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
                                          #HLTPaths = cms.vstring('HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v*','HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg_v*','HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v*'),
                                          HLTPaths = cms.vstring('*'),
                                          andOr = cms.bool(True),
                                          throw = cms.bool(False)
                                          )

process.p = cms.Path(
process.siStripMatchedRecHits*process.siPixelRecHits*process.MeasurementTrackerEvent*process.TrackRefitter*
#  process.hltFilter*
#  process.patDefaultSequence*
process.fevt
)
