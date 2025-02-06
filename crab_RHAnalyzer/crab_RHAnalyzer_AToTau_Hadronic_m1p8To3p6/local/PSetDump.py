import FWCore.ParameterSet.Config as cms
from HeterogeneousCore.CUDACore.SwitchProducerCUDA import SwitchProducerCUDA

process = cms.Process("FEVTAnalyzer")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(),
    skipEvents = cms.untracked.uint32(0)
)
process.ChargeSignificanceTrajectoryFilter_block = cms.PSet(
    ComponentType = cms.string('ChargeSignificanceTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0)
)

process.CkfBaseTrajectoryFilter_block = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    highEtaSwitch = cms.double(5.0),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(999),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutNone')
    ),
    minHitsAtHighEta = cms.int32(5),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.9),
    minimumNumberOfHits = cms.int32(5),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.CkfTrajectoryBuilder = cms.PSet(
    ComponentType = cms.string('CkfTrajectoryBuilder'),
    TTRHBuilder = cms.string('WithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    estimator = cms.string('Chi2'),
    intermediateCleaning = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterial'),
    propagatorOpposite = cms.string('PropagatorWithMaterialOpposite'),
    seedAs5DHit = cms.bool(False),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('CkfBaseTrajectoryFilter_block')
    ),
    updator = cms.string('KFUpdator')
)

process.CkfTrajectoryBuilderBeamHalo = cms.PSet(
    ComponentType = cms.string('CkfTrajectoryBuilder'),
    TTRHBuilder = cms.string('WithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    estimator = cms.string('Chi2'),
    intermediateCleaning = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(5),
    propagatorAlong = cms.string('BeamHaloPropagatorAlong'),
    propagatorOpposite = cms.string('BeamHaloPropagatorOpposite'),
    seedAs5DHit = cms.bool(False),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('ckfTrajectoryFilterBeamHaloMuon')
    ),
    updator = cms.string('KFUpdator')
)

process.ClusterShapeTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('ClusterShapeTrajectoryFilter'),
    cacheSrc = cms.InputTag("siPixelClusterShapeCache")
)

process.CommonClusterCheckPSet = cms.PSet(
    ClusterCollectionLabel = cms.InputTag("siStripClusters"),
    MaxNumberOfCosmicClusters = cms.uint32(400000),
    MaxNumberOfPixelClusters = cms.uint32(40000),
    PixelClusterCollectionLabel = cms.InputTag("siPixelClusters"),
    cut = cms.string('strip < 400000 && pixel < 40000 && (strip < 50000 + 10*pixel) && (pixel < 5000 + 0.1*strip)'),
    doClusterCheck = cms.bool(False)
)

process.CompositeTrajectoryFilter_block = cms.PSet(
    ComponentType = cms.string('CompositeTrajectoryFilter'),
    filters = cms.VPSet()
)

process.CondDB = cms.PSet(
    DBParameters = cms.PSet(
        authenticationPath = cms.untracked.string(''),
        authenticationSystem = cms.untracked.int32(0),
        connectionTimeout = cms.untracked.int32(0),
        messageLevel = cms.untracked.int32(0),
        security = cms.untracked.string('')
    ),
    connect = cms.string('')
)

process.CosmicSeedCreator = cms.PSet(
    ComponentName = cms.string('CosmicSeedCreator'),
    MinOneOverPtError = cms.double(1.0),
    OriginTransverseErrorMultiplier = cms.double(1.0),
    SeedMomentumForBOFF = cms.double(5.0),
    TTRHBuilder = cms.string('WithTrackAngle'),
    forceKinematicWithRegionDirection = cms.bool(False),
    magneticField = cms.string('ParabolicMf'),
    maxseeds = cms.int32(10000),
    propagator = cms.string('PropagatorWithMaterialParabolicMf')
)

process.DefaultAlgorithms = cms.PSet(
    APVInspectMode = cms.string('BaselineFollower'),
    APVRestoreMode = cms.string('BaselineFollower'),
    ApplyBaselineCleaner = cms.bool(True),
    ApplyBaselineRejection = cms.bool(True),
    CleaningSequence = cms.uint32(1),
    CommonModeNoiseSubtractionMode = cms.string('IteratedMedian'),
    CutToAvoidSignal = cms.double(2.0),
    DeltaCMThreshold = cms.uint32(20),
    Deviation = cms.uint32(25),
    ForceNoRestore = cms.bool(False),
    Fraction = cms.double(0.2),
    Iterations = cms.int32(3),
    MeanCM = cms.int32(0),
    PedestalSubtractionFedMode = cms.bool(False),
    SiStripFedZeroSuppressionMode = cms.uint32(4),
    TruncateInSuppressor = cms.bool(True),
    Use10bitsTruncation = cms.bool(False),
    consecThreshold = cms.uint32(5),
    discontinuityThreshold = cms.int32(12),
    distortionThreshold = cms.uint32(20),
    doAPVRestore = cms.bool(True),
    filteredBaselineDerivativeSumSquare = cms.double(30),
    filteredBaselineMax = cms.double(6),
    hitStripThreshold = cms.uint32(40),
    lastGradient = cms.int32(10),
    minStripsToFit = cms.uint32(4),
    nSaturatedStrip = cms.uint32(2),
    nSigmaNoiseDerTh = cms.uint32(4),
    nSmooth = cms.uint32(9),
    restoreThreshold = cms.double(0.5),
    sizeWindow = cms.int32(1),
    slopeX = cms.int32(3),
    slopeY = cms.int32(4),
    useCMMeanMap = cms.bool(False),
    useRealMeanCM = cms.bool(False),
    widthCluster = cms.int32(64)
)

process.DefaultClusterizer = cms.PSet(
    Algorithm = cms.string('ThreeThresholdAlgorithm'),
    ChannelThreshold = cms.double(2.0),
    ClusterThreshold = cms.double(5.0),
    ConditionsLabel = cms.string(''),
    MaxAdjacentBad = cms.uint32(0),
    MaxSequentialBad = cms.uint32(1),
    MaxSequentialHoles = cms.uint32(0),
    RemoveApvShots = cms.bool(True),
    SeedThreshold = cms.double(3.0),
    clusterChargeCut = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutNone')
    )
)

process.GroupedCkfTrajectoryBuilder = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    TTRHBuilder = cms.string('WithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('Chi2'),
    foundHitBonus = cms.double(10.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('CkfBaseTrajectoryFilter_block')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(5),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.0),
    minNrOfHitsForRebuild = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterial'),
    propagatorOpposite = cms.string('PropagatorWithMaterialOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    seedAs5DHit = cms.bool(False),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('CkfBaseTrajectoryFilter_block')
    ),
    updator = cms.string('KFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.GroupedCkfTrajectoryBuilderCDC = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    TTRHBuilder = cms.string('WithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('Chi2MeasurementEstimatorForCDC'),
    foundHitBonus = cms.double(10.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('CkfBaseTrajectoryFilter_block')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(3),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.0),
    minNrOfHitsForRebuild = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterial'),
    propagatorOpposite = cms.string('PropagatorWithMaterialOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    seedAs5DHit = cms.bool(False),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('ckfBaseTrajectoryFilterCDC')
    ),
    updator = cms.string('KFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.GroupedCkfTrajectoryBuilderP5 = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    TTRHBuilder = cms.string('WithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('Chi2MeasurementEstimatorForP5'),
    foundHitBonus = cms.double(10.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('CkfBaseTrajectoryFilter_block')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(1),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.0),
    minNrOfHitsForRebuild = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterial'),
    propagatorOpposite = cms.string('PropagatorWithMaterialOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    seedAs5DHit = cms.bool(False),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('ckfBaseTrajectoryFilterP5')
    ),
    updator = cms.string('KFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.GroupedCkfTrajectoryBuilderP5Bottom = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    TTRHBuilder = cms.string('WithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('Chi2MeasurementEstimatorForP5'),
    foundHitBonus = cms.double(10.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('CkfBaseTrajectoryFilter_block')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(1),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.0),
    minNrOfHitsForRebuild = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterial'),
    propagatorOpposite = cms.string('PropagatorWithMaterialOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    seedAs5DHit = cms.bool(False),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('ckfBaseTrajectoryFilterP5')
    ),
    updator = cms.string('KFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.GroupedCkfTrajectoryBuilderP5Top = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    TTRHBuilder = cms.string('WithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('Chi2MeasurementEstimatorForP5'),
    foundHitBonus = cms.double(10.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('CkfBaseTrajectoryFilter_block')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(1),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.0),
    minNrOfHitsForRebuild = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterial'),
    propagatorOpposite = cms.string('PropagatorWithMaterialOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    seedAs5DHit = cms.bool(False),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('ckfBaseTrajectoryFilterP5')
    ),
    updator = cms.string('KFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.HFRecalParameterBlock = cms.PSet(
    HFdepthOneParameterA = cms.vdouble(
        0.004123, 0.00602, 0.008201, 0.010489, 0.013379,
        0.016997, 0.021464, 0.027371, 0.034195, 0.044807,
        0.058939, 0.125497
    ),
    HFdepthOneParameterB = cms.vdouble(
        -4e-06, -2e-06, 0.0, 4e-06, 1.5e-05,
        2.6e-05, 6.3e-05, 8.4e-05, 0.00016, 0.000107,
        0.000425, 0.000209
    ),
    HFdepthTwoParameterA = cms.vdouble(
        0.002861, 0.004168, 0.0064, 0.008388, 0.011601,
        0.014425, 0.018633, 0.023232, 0.028274, 0.035447,
        0.051579, 0.086593
    ),
    HFdepthTwoParameterB = cms.vdouble(
        -2e-06, -0.0, -7e-06, -6e-06, -2e-06,
        1e-06, 1.9e-05, 3.1e-05, 6.7e-05, 1.2e-05,
        0.000157, -3e-06
    )
)

process.METSignificanceParams = cms.PSet(
    dRMatch = cms.double(0.4),
    jetThreshold = cms.double(15),
    jeta = cms.vdouble(0.8, 1.3, 1.9, 2.5),
    jpar = cms.vdouble(1.39, 1.26, 1.21, 1.23, 1.28),
    pjpar = cms.vdouble(-0.2586, 0.6173)
)

process.MaxCCCLostHitsTrajectoryFilter_block = cms.PSet(
    ComponentType = cms.string('MaxCCCLostHitsTrajectoryFilter'),
    maxCCCLostHits = cms.int32(3),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutLoose')
    )
)

process.MaxConsecLostHitsTrajectoryFilter_block = cms.PSet(
    ComponentType = cms.string('MaxConsecLostHitsTrajectoryFilter'),
    maxConsecLostHits = cms.int32(1)
)

process.MaxHitsTrajectoryFilter_block = cms.PSet(
    ComponentType = cms.string('MaxHitsTrajectoryFilter'),
    maxNumberOfHits = cms.int32(100)
)

process.MaxLostHitsTrajectoryFilter_block = cms.PSet(
    ComponentType = cms.string('MaxLostHitsTrajectoryFilter'),
    maxLostHits = cms.int32(2)
)

process.MinHitsTrajectoryFilter_block = cms.PSet(
    ComponentType = cms.string('MinHitsTrajectoryFilter'),
    highEtaSwitch = cms.double(5.0),
    minHitsAtHighEta = cms.int32(5),
    minimumNumberOfHits = cms.int32(5)
)

process.MinPtTrajectoryFilter_block = cms.PSet(
    ComponentType = cms.string('MinPtTrajectoryFilter'),
    minHitsMinPt = cms.int32(3),
    minPt = cms.double(1.0),
    nSigmaMinPt = cms.double(5.0)
)

process.PFRecoTauPFJetInputs = cms.PSet(
    inputJetCollection = cms.InputTag("ak4PFJets"),
    isolationConeSize = cms.double(0.5),
    jetConeSize = cms.double(0.5),
    maxJetAbsEta = cms.double(2.7),
    minJetPt = cms.double(14.0)
)

process.PFTauQualityCuts = cms.PSet(
    isolationQualityCuts = cms.PSet(
        maxDeltaZ = cms.double(0.2),
        maxDeltaZToLeadTrack = cms.double(-1.0),
        maxTrackChi2 = cms.double(100.0),
        maxTransverseImpactParameter = cms.double(0.03),
        minGammaEt = cms.double(1.5),
        minTrackHits = cms.uint32(8),
        minTrackPixelHits = cms.uint32(0),
        minTrackPt = cms.double(1.0),
        minTrackVertexWeight = cms.double(-1.0)
    ),
    leadingTrkOrPFCandOption = cms.string('leadPFCand'),
    primaryVertexSrc = cms.InputTag("offlinePrimaryVertices"),
    pvFindingAlgo = cms.string('closestInDeltaZ'),
    recoverLeadingTrk = cms.bool(False),
    signalQualityCuts = cms.PSet(
        maxDeltaZ = cms.double(0.4),
        maxDeltaZToLeadTrack = cms.double(-1.0),
        maxTrackChi2 = cms.double(100.0),
        maxTransverseImpactParameter = cms.double(0.1),
        minGammaEt = cms.double(1.0),
        minNeutralHadronEt = cms.double(30.0),
        minTrackHits = cms.uint32(3),
        minTrackPixelHits = cms.uint32(0),
        minTrackPt = cms.double(0.5),
        minTrackVertexWeight = cms.double(-1.0)
    ),
    vertexTrackFiltering = cms.bool(False),
    vxAssocQualityCuts = cms.PSet(
        maxTrackChi2 = cms.double(100.0),
        maxTransverseImpactParameter = cms.double(0.1),
        minGammaEt = cms.double(1.0),
        minTrackHits = cms.uint32(3),
        minTrackPixelHits = cms.uint32(0),
        minTrackPt = cms.double(0.5),
        minTrackVertexWeight = cms.double(-1.0)
    )
)

process.PixelTripletHLTGenerator = cms.PSet(
    ComponentName = cms.string('PixelTripletHLTGenerator'),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('none')
    ),
    extraHitRPhitolerance = cms.double(0.032),
    extraHitRZtolerance = cms.double(0.037),
    maxElement = cms.uint32(100000),
    phiPreFiltering = cms.double(0.3),
    useBending = cms.bool(True),
    useFixedPreFiltering = cms.bool(False),
    useMultScattering = cms.bool(True)
)

process.PixelTripletHLTGeneratorWithFilter = cms.PSet(
    ComponentName = cms.string('PixelTripletHLTGenerator'),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('LowPtClusterShapeSeedComparitor'),
        clusterShapeCacheSrc = cms.InputTag("siPixelClusterShapeCache"),
        clusterShapeHitFilter = cms.string('ClusterShapeHitFilter')
    ),
    extraHitRPhitolerance = cms.double(0.032),
    extraHitRZtolerance = cms.double(0.037),
    maxElement = cms.uint32(100000),
    phiPreFiltering = cms.double(0.3),
    useBending = cms.bool(True),
    useFixedPreFiltering = cms.bool(False),
    useMultScattering = cms.bool(True)
)

process.RegionPSetBlock = cms.PSet(
    RegionPSet = cms.PSet(
        originHalfLength = cms.double(21.2),
        originRadius = cms.double(0.2),
        originXPos = cms.double(0.0),
        originYPos = cms.double(0.0),
        originZPos = cms.double(0.0),
        precise = cms.bool(True),
        ptMin = cms.double(0.9),
        useMultipleScattering = cms.bool(False)
    )
)

process.SiStripClusterChargeCutLoose = cms.PSet(
    value = cms.double(1620.0)
)

process.SiStripClusterChargeCutNone = cms.PSet(
    value = cms.double(-1.0)
)

process.SiStripClusterChargeCutTight = cms.PSet(
    value = cms.double(1945.0)
)

process.SiStripClusterChargeCutTiny = cms.PSet(
    value = cms.double(800.0)
)

process.ThresholdPtTrajectoryFilter_block = cms.PSet(
    ComponentType = cms.string('ThresholdPtTrajectoryFilter'),
    minHitsThresholdPt = cms.int32(3),
    nSigmaThresholdPt = cms.double(5.0),
    thresholdPt = cms.double(10.0)
)

process.ckfBaseInOutTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    highEtaSwitch = cms.double(5.0),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(999),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutNone')
    ),
    minHitsAtHighEta = cms.int32(5),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.9),
    minimumNumberOfHits = cms.int32(5),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.ckfBaseTrajectoryFilterCDC = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    highEtaSwitch = cms.double(5.0),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(10),
    maxLostHits = cms.int32(10),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutNone')
    ),
    minHitsAtHighEta = cms.int32(5),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.5),
    minimumNumberOfHits = cms.int32(5),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.ckfBaseTrajectoryFilterP5 = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    highEtaSwitch = cms.double(5.0),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(3),
    maxLostHits = cms.int32(4),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutNone')
    ),
    minHitsAtHighEta = cms.int32(5),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.5),
    minimumNumberOfHits = cms.int32(5),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.ckfTrajectoryFilterBeamHaloMuon = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    highEtaSwitch = cms.double(5.0),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(2),
    maxLostHits = cms.int32(3),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutNone')
    ),
    minHitsAtHighEta = cms.int32(5),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.1),
    minimumNumberOfHits = cms.int32(4),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.conv2CkfTrajectoryBuilder = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    TTRHBuilder = cms.string('WithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('Chi2'),
    foundHitBonus = cms.double(10.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('CkfBaseTrajectoryFilter_block')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(2),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.0),
    minNrOfHitsForRebuild = cms.int32(3),
    propagatorAlong = cms.string('PropagatorWithMaterial'),
    propagatorOpposite = cms.string('PropagatorWithMaterialOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    seedAs5DHit = cms.bool(False),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('conv2CkfTrajectoryFilter')
    ),
    updator = cms.string('KFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.conv2CkfTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    highEtaSwitch = cms.double(5.0),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(1),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutNone')
    ),
    minHitsAtHighEta = cms.int32(5),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.1),
    minimumNumberOfHits = cms.int32(3),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.convCkfTrajectoryBuilder = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    TTRHBuilder = cms.string('WithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('convStepChi2Est'),
    foundHitBonus = cms.double(10.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('CkfBaseTrajectoryFilter_block')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(1),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.0),
    minNrOfHitsForRebuild = cms.int32(3),
    propagatorAlong = cms.string('PropagatorWithMaterial'),
    propagatorOpposite = cms.string('PropagatorWithMaterialOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    seedAs5DHit = cms.bool(False),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('convCkfTrajectoryFilter')
    ),
    updator = cms.string('KFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.convCkfTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    highEtaSwitch = cms.double(5.0),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(1),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutNone')
    ),
    minHitsAtHighEta = cms.int32(5),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.1),
    minimumNumberOfHits = cms.int32(3),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.decayMode_1Prong0Pi0 = cms.PSet(
    applyBendCorrection = cms.PSet(
        eta = cms.bool(True),
        mass = cms.bool(True),
        phi = cms.bool(True)
    ),
    assumeStripMass = cms.double(-1.0),
    maxMass = cms.string('1.'),
    maxPi0Mass = cms.double(1000000000.0),
    minMass = cms.double(-1000.0),
    minPi0Mass = cms.double(-1000.0),
    nCharged = cms.uint32(1),
    nChargedPFCandsMin = cms.uint32(1),
    nPiZeros = cms.uint32(0),
    nTracksMin = cms.uint32(1)
)

process.decayMode_1Prong1Pi0 = cms.PSet(
    applyBendCorrection = cms.PSet(
        eta = cms.bool(True),
        mass = cms.bool(True),
        phi = cms.bool(True)
    ),
    assumeStripMass = cms.double(0.1349),
    maxMass = cms.string('max(1.3, min(1.3*sqrt(pt/100.), 4.2))'),
    maxPi0Mass = cms.double(1000000000.0),
    minMass = cms.double(0.3),
    minPi0Mass = cms.double(-1000.0),
    nCharged = cms.uint32(1),
    nChargedPFCandsMin = cms.uint32(1),
    nPiZeros = cms.uint32(1),
    nTracksMin = cms.uint32(1)
)

process.decayMode_1Prong2Pi0 = cms.PSet(
    applyBendCorrection = cms.PSet(
        eta = cms.bool(True),
        mass = cms.bool(True),
        phi = cms.bool(True)
    ),
    assumeStripMass = cms.double(0.0),
    maxMass = cms.string('max(1.2, min(1.2*sqrt(pt/100.), 4.0))'),
    maxPi0Mass = cms.double(0.2),
    minMass = cms.double(0.4),
    minPi0Mass = cms.double(0.05),
    nCharged = cms.uint32(1),
    nChargedPFCandsMin = cms.uint32(1),
    nPiZeros = cms.uint32(2),
    nTracksMin = cms.uint32(1)
)

process.decayMode_2Prong0Pi0 = cms.PSet(
    applyBendCorrection = cms.PSet(
        eta = cms.bool(False),
        mass = cms.bool(False),
        phi = cms.bool(False)
    ),
    assumeStripMass = cms.double(-1.0),
    maxMass = cms.string('1.2'),
    maxPi0Mass = cms.double(1000000000.0),
    minMass = cms.double(0.0),
    minPi0Mass = cms.double(-1000.0),
    nCharged = cms.uint32(2),
    nChargedPFCandsMin = cms.uint32(1),
    nPiZeros = cms.uint32(0),
    nTracksMin = cms.uint32(2)
)

process.decayMode_2Prong1Pi0 = cms.PSet(
    applyBendCorrection = cms.PSet(
        eta = cms.bool(False),
        mass = cms.bool(False),
        phi = cms.bool(False)
    ),
    assumeStripMass = cms.double(-1.0),
    maxMass = cms.string('max(1.2, min(1.2*sqrt(pt/100.), 4.0))'),
    maxPi0Mass = cms.double(1000000000.0),
    minMass = cms.double(0.0),
    minPi0Mass = cms.double(-1000.0),
    nCharged = cms.uint32(2),
    nChargedPFCandsMin = cms.uint32(1),
    nPiZeros = cms.uint32(1),
    nTracksMin = cms.uint32(2)
)

process.decayMode_3Prong0Pi0 = cms.PSet(
    applyBendCorrection = cms.PSet(
        eta = cms.bool(False),
        mass = cms.bool(False),
        phi = cms.bool(False)
    ),
    assumeStripMass = cms.double(-1.0),
    maxMass = cms.string('1.5'),
    maxPi0Mass = cms.double(1000000000.0),
    minMass = cms.double(0.8),
    minPi0Mass = cms.double(-1000.0),
    nCharged = cms.uint32(3),
    nChargedPFCandsMin = cms.uint32(1),
    nPiZeros = cms.uint32(0),
    nTracksMin = cms.uint32(2)
)

process.decayMode_3Prong1Pi0 = cms.PSet(
    applyBendCorrection = cms.PSet(
        eta = cms.bool(False),
        mass = cms.bool(False),
        phi = cms.bool(False)
    ),
    assumeStripMass = cms.double(-1.0),
    maxMass = cms.string('1.6'),
    maxPi0Mass = cms.double(1000000000.0),
    minMass = cms.double(0.9),
    minPi0Mass = cms.double(-1000.0),
    nCharged = cms.uint32(3),
    nChargedPFCandsMin = cms.uint32(1),
    nPiZeros = cms.uint32(1),
    nTracksMin = cms.uint32(2)
)

process.detachedQuadStepTrajectoryBuilder = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    TTRHBuilder = cms.string('WithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('detachedQuadStepChi2Est'),
    foundHitBonus = cms.double(10.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('CkfBaseTrajectoryFilter_block')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(3),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.0),
    minNrOfHitsForRebuild = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterial'),
    propagatorOpposite = cms.string('PropagatorWithMaterialOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    seedAs5DHit = cms.bool(False),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('detachedQuadStepTrajectoryFilter')
    ),
    updator = cms.string('KFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.detachedQuadStepTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('CompositeTrajectoryFilter'),
    filters = cms.VPSet(cms.PSet(
        refToPSet_ = cms.string('detachedQuadStepTrajectoryFilterBase')
    ))
)

process.detachedQuadStepTrajectoryFilterBase = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    highEtaSwitch = cms.double(5.0),
    maxCCCLostHits = cms.int32(0),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(999),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutLoose')
    ),
    minHitsAtHighEta = cms.int32(5),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.075),
    minimumNumberOfHits = cms.int32(3),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.detachedTripletStepTrajectoryBuilder = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    TTRHBuilder = cms.string('WithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('detachedTripletStepChi2Est'),
    foundHitBonus = cms.double(10.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('CkfBaseTrajectoryFilter_block')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(3),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.0),
    minNrOfHitsForRebuild = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterial'),
    propagatorOpposite = cms.string('PropagatorWithMaterialOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    seedAs5DHit = cms.bool(False),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('detachedTripletStepTrajectoryFilter')
    ),
    updator = cms.string('KFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.detachedTripletStepTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('CompositeTrajectoryFilter'),
    filters = cms.VPSet(cms.PSet(
        refToPSet_ = cms.string('detachedTripletStepTrajectoryFilterBase')
    ))
)

process.detachedTripletStepTrajectoryFilterBase = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    highEtaSwitch = cms.double(5.0),
    maxCCCLostHits = cms.int32(0),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(999),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutLoose')
    ),
    minHitsAtHighEta = cms.int32(5),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.075),
    minimumNumberOfHits = cms.int32(3),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.detachedTripletStepTrajectoryFilterShape = cms.PSet(
    ComponentType = cms.string('StripSubClusterShapeTrajectoryFilter'),
    layerMask = cms.PSet(
        TEC = cms.bool(False),
        TIB = cms.vuint32(1, 2),
        TID = cms.vuint32(1, 2),
        TOB = cms.bool(False)
    ),
    maxNSat = cms.uint32(3),
    maxTrimmedSizeDiffNeg = cms.double(1.0),
    maxTrimmedSizeDiffPos = cms.double(0.7),
    seedCutMIPs = cms.double(0.35),
    seedCutSN = cms.double(7.0),
    subclusterCutMIPs = cms.double(0.45),
    subclusterCutSN = cms.double(12.0),
    subclusterWindow = cms.double(0.7),
    trimMaxADC = cms.double(30.0),
    trimMaxFracNeigh = cms.double(0.25),
    trimMaxFracTotal = cms.double(0.15)
)

process.displacedGeneralStepTrajectoryBuilder = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    TTRHBuilder = cms.string('WithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('displacedGeneralStepChi2Est'),
    foundHitBonus = cms.double(10.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('displacedGeneralStepTrajectoryFilterInOut')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(2),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.0),
    minNrOfHitsForRebuild = cms.int32(4),
    propagatorAlong = cms.string('PropagatorWithMaterial'),
    propagatorOpposite = cms.string('PropagatorWithMaterialOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    seedAs5DHit = cms.bool(False),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('displacedGeneralStepTrajectoryFilter')
    ),
    updator = cms.string('KFUpdator'),
    useSameTrajFilter = cms.bool(False)
)

process.displacedGeneralStepTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    highEtaSwitch = cms.double(5.0),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(1),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutNone')
    ),
    minHitsAtHighEta = cms.int32(5),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(1),
    minimumNumberOfHits = cms.int32(4),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(1),
    strictSeedExtension = cms.bool(False)
)

process.displacedGeneralStepTrajectoryFilterInOut = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    highEtaSwitch = cms.double(5.0),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(1),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutNone')
    ),
    minHitsAtHighEta = cms.int32(5),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(1),
    minimumNumberOfHits = cms.int32(4),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(1),
    strictSeedExtension = cms.bool(False)
)

process.highPtTripletStepTrajectoryBuilder = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    TTRHBuilder = cms.string('WithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('highPtTripletStepChi2Est'),
    foundHitBonus = cms.double(10.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('CkfBaseTrajectoryFilter_block')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(3),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.0),
    minNrOfHitsForRebuild = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterial'),
    propagatorOpposite = cms.string('PropagatorWithMaterialOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    seedAs5DHit = cms.bool(False),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('highPtTripletStepTrajectoryFilter')
    ),
    updator = cms.string('KFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.highPtTripletStepTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('CompositeTrajectoryFilter'),
    filters = cms.VPSet(cms.PSet(
        refToPSet_ = cms.string('highPtTripletStepTrajectoryFilterBase')
    ))
)

process.highPtTripletStepTrajectoryFilterBase = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    highEtaSwitch = cms.double(5.0),
    maxCCCLostHits = cms.int32(0),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(999),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutLoose')
    ),
    minHitsAtHighEta = cms.int32(5),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.2),
    minimumNumberOfHits = cms.int32(3),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.highPtTripletStepTrajectoryFilterInOut = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    highEtaSwitch = cms.double(5.0),
    maxCCCLostHits = cms.int32(0),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(999),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutLoose')
    ),
    minHitsAtHighEta = cms.int32(5),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.4),
    minimumNumberOfHits = cms.int32(4),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(1),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.initialStepTrajectoryBuilder = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    TTRHBuilder = cms.string('WithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('initialStepChi2Est'),
    foundHitBonus = cms.double(10.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('CkfBaseTrajectoryFilter_block')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(True),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(3),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.0),
    minNrOfHitsForRebuild = cms.int32(1),
    propagatorAlong = cms.string('PropagatorWithMaterial'),
    propagatorOpposite = cms.string('PropagatorWithMaterialOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    seedAs5DHit = cms.bool(False),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('initialStepTrajectoryFilter')
    ),
    updator = cms.string('KFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.initialStepTrajectoryBuilderPreSplitting = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    TTRHBuilder = cms.string('WithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('initialStepChi2Est'),
    foundHitBonus = cms.double(10.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('CkfBaseTrajectoryFilter_block')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(3),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.0),
    minNrOfHitsForRebuild = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterial'),
    propagatorOpposite = cms.string('PropagatorWithMaterialOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    seedAs5DHit = cms.bool(False),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('initialStepTrajectoryFilterPreSplitting')
    ),
    updator = cms.string('KFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.initialStepTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('CompositeTrajectoryFilter'),
    filters = cms.VPSet(cms.PSet(
        refToPSet_ = cms.string('initialStepTrajectoryFilterBase')
    ))
)

process.initialStepTrajectoryFilterBase = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    highEtaSwitch = cms.double(5.0),
    maxCCCLostHits = cms.int32(0),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(999),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutLoose')
    ),
    minHitsAtHighEta = cms.int32(5),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.2),
    minimumNumberOfHits = cms.int32(3),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.initialStepTrajectoryFilterBasePreSplitting = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    highEtaSwitch = cms.double(5.0),
    maxCCCLostHits = cms.int32(0),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(999),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutLoose')
    ),
    minHitsAtHighEta = cms.int32(5),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.2),
    minimumNumberOfHits = cms.int32(4),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.initialStepTrajectoryFilterInOut = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    highEtaSwitch = cms.double(5.0),
    maxCCCLostHits = cms.int32(0),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(999),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutLoose')
    ),
    minHitsAtHighEta = cms.int32(5),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.2),
    minimumNumberOfHits = cms.int32(4),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(True),
    seedExtension = cms.int32(1),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(True)
)

process.initialStepTrajectoryFilterPreSplitting = cms.PSet(
    ComponentType = cms.string('CompositeTrajectoryFilter'),
    filters = cms.VPSet(
        cms.PSet(
            refToPSet_ = cms.string('initialStepTrajectoryFilterBasePreSplitting')
        ),
        cms.PSet(
            refToPSet_ = cms.string('initialStepTrajectoryFilterShapePreSplitting')
        )
    )
)

process.initialStepTrajectoryFilterShape = cms.PSet(
    ComponentType = cms.string('StripSubClusterShapeTrajectoryFilter'),
    layerMask = cms.PSet(
        TEC = cms.bool(False),
        TIB = cms.vuint32(1, 2),
        TID = cms.vuint32(1, 2),
        TOB = cms.bool(False)
    ),
    maxNSat = cms.uint32(3),
    maxTrimmedSizeDiffNeg = cms.double(1.0),
    maxTrimmedSizeDiffPos = cms.double(0.7),
    seedCutMIPs = cms.double(0.35),
    seedCutSN = cms.double(7.0),
    subclusterCutMIPs = cms.double(0.45),
    subclusterCutSN = cms.double(12.0),
    subclusterWindow = cms.double(0.7),
    trimMaxADC = cms.double(30.0),
    trimMaxFracNeigh = cms.double(0.25),
    trimMaxFracTotal = cms.double(0.15)
)

process.initialStepTrajectoryFilterShapePreSplitting = cms.PSet(
    ComponentType = cms.string('StripSubClusterShapeTrajectoryFilter'),
    layerMask = cms.PSet(
        TEC = cms.bool(False),
        TIB = cms.vuint32(1, 2),
        TID = cms.vuint32(1, 2),
        TOB = cms.bool(False)
    ),
    maxNSat = cms.uint32(3),
    maxTrimmedSizeDiffNeg = cms.double(1.0),
    maxTrimmedSizeDiffPos = cms.double(0.7),
    seedCutMIPs = cms.double(0.35),
    seedCutSN = cms.double(7.0),
    subclusterCutMIPs = cms.double(0.45),
    subclusterCutSN = cms.double(12.0),
    subclusterWindow = cms.double(0.7),
    trimMaxADC = cms.double(30.0),
    trimMaxFracNeigh = cms.double(0.25),
    trimMaxFracTotal = cms.double(0.15)
)

process.jetCoreRegionalStepBarrelTrajectoryBuilder = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    TTRHBuilder = cms.string('WithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('jetCoreRegionalStepChi2Est'),
    foundHitBonus = cms.double(10.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('CkfBaseTrajectoryFilter_block')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(True),
    lockHits = cms.bool(False),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(50),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.0),
    minNrOfHitsForRebuild = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterial'),
    propagatorOpposite = cms.string('PropagatorWithMaterialOpposite'),
    requireSeedHitsInRebuild = cms.bool(False),
    seedAs5DHit = cms.bool(False),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('jetCoreRegionalStepBarrelTrajectoryFilter')
    ),
    updator = cms.string('KFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.jetCoreRegionalStepBarrelTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    highEtaSwitch = cms.double(5.0),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(2),
    maxLostHits = cms.int32(999),
    maxLostHitsFraction = cms.double(1.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutNone')
    ),
    minHitsAtHighEta = cms.int32(5),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.9),
    minimumNumberOfHits = cms.int32(2),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.jetCoreRegionalStepEndcapTrajectoryBuilder = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    TTRHBuilder = cms.string('WithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('jetCoreRegionalStepChi2Est'),
    foundHitBonus = cms.double(10.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('CkfBaseTrajectoryFilter_block')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(50),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.0),
    minNrOfHitsForRebuild = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterial'),
    propagatorOpposite = cms.string('PropagatorWithMaterialOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    seedAs5DHit = cms.bool(False),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('jetCoreRegionalStepEndcapTrajectoryFilter')
    ),
    updator = cms.string('KFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.jetCoreRegionalStepEndcapTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    highEtaSwitch = cms.double(5.0),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(999),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutNone')
    ),
    minHitsAtHighEta = cms.int32(5),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.1),
    minimumNumberOfHits = cms.int32(4),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.jetCoreRegionalStepTrajectoryBuilder = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    TTRHBuilder = cms.string('WithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('jetCoreRegionalStepChi2Est'),
    foundHitBonus = cms.double(10.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('CkfBaseTrajectoryFilter_block')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(50),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.0),
    minNrOfHitsForRebuild = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterial'),
    propagatorOpposite = cms.string('PropagatorWithMaterialOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    seedAs5DHit = cms.bool(False),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('jetCoreRegionalStepTrajectoryFilter')
    ),
    updator = cms.string('KFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.jetCoreRegionalStepTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    highEtaSwitch = cms.double(5.0),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(999),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutNone')
    ),
    minHitsAtHighEta = cms.int32(5),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.1),
    minimumNumberOfHits = cms.int32(4),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.layerInfo = cms.PSet(
    FPix = cms.PSet(
        HitProducer = cms.string('siPixelRecHits'),
        TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4PixelPairs')
    ),
    MTEC = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        maxRing = cms.int32(7),
        minRing = cms.int32(5),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
        useRingSlector = cms.bool(False)
    ),
    TEC = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        maxRing = cms.int32(2),
        minRing = cms.int32(1),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
        useRingSlector = cms.bool(True)
    ),
    TID = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
        useRingSlector = cms.bool(False)
    )
)

process.leadTrackFinding = cms.PSet(
    Producer = cms.InputTag("pfRecoTauDiscriminationByLeadingTrackFinding"),
    cut = cms.double(0.5)
)

process.lowPtQuadStepTrajectoryBuilder = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    TTRHBuilder = cms.string('WithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('lowPtQuadStepChi2Est'),
    foundHitBonus = cms.double(10.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('CkfBaseTrajectoryFilter_block')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(4),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.0),
    minNrOfHitsForRebuild = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterial'),
    propagatorOpposite = cms.string('PropagatorWithMaterialOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    seedAs5DHit = cms.bool(False),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('lowPtQuadStepTrajectoryFilter')
    ),
    updator = cms.string('KFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.lowPtQuadStepTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('CompositeTrajectoryFilter'),
    filters = cms.VPSet(cms.PSet(
        refToPSet_ = cms.string('lowPtQuadStepTrajectoryFilterBase')
    ))
)

process.lowPtQuadStepTrajectoryFilterBase = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    highEtaSwitch = cms.double(5.0),
    maxCCCLostHits = cms.int32(0),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(999),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutLoose')
    ),
    minHitsAtHighEta = cms.int32(5),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.075),
    minimumNumberOfHits = cms.int32(3),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.lowPtTripletStepStandardTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    highEtaSwitch = cms.double(5.0),
    maxCCCLostHits = cms.int32(0),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(999),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutLoose')
    ),
    minHitsAtHighEta = cms.int32(5),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.075),
    minimumNumberOfHits = cms.int32(3),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.lowPtTripletStepTrajectoryBuilder = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    TTRHBuilder = cms.string('WithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('lowPtTripletStepChi2Est'),
    foundHitBonus = cms.double(10.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('CkfBaseTrajectoryFilter_block')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(4),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.0),
    minNrOfHitsForRebuild = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterial'),
    propagatorOpposite = cms.string('PropagatorWithMaterialOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    seedAs5DHit = cms.bool(False),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('lowPtTripletStepTrajectoryFilter')
    ),
    updator = cms.string('KFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.lowPtTripletStepTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('CompositeTrajectoryFilter'),
    filters = cms.VPSet(cms.PSet(
        refToPSet_ = cms.string('lowPtTripletStepStandardTrajectoryFilter')
    ))
)

process.lowPtTripletStepTrajectoryFilterInOut = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    highEtaSwitch = cms.double(5.0),
    maxCCCLostHits = cms.int32(0),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(999),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutLoose')
    ),
    minHitsAtHighEta = cms.int32(5),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.075),
    minimumNumberOfHits = cms.int32(4),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(1),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

process.maxLuminosityBlocks = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.mixedTripletStepTrajectoryBuilder = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    TTRHBuilder = cms.string('WithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('mixedTripletStepChi2Est'),
    foundHitBonus = cms.double(10.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('CkfBaseTrajectoryFilter_block')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(2),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.0),
    minNrOfHitsForRebuild = cms.int32(5),
    propagatorAlong = cms.string('mixedTripletStepPropagator'),
    propagatorOpposite = cms.string('mixedTripletStepPropagatorOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    seedAs5DHit = cms.bool(False),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('mixedTripletStepTrajectoryFilter')
    ),
    updator = cms.string('KFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.mixedTripletStepTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(1.4),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    highEtaSwitch = cms.double(5.0),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(999),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutNone')
    ),
    minHitsAtHighEta = cms.int32(5),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.1),
    minimumNumberOfHits = cms.int32(3),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.muonSeededTrajectoryBuilderForInOut = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    TTRHBuilder = cms.string('WithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('muonSeededMeasurementEstimatorForInOut'),
    foundHitBonus = cms.double(1000.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('muonSeededTrajectoryFilterForInOut')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(True),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(1.0),
    maxCand = cms.int32(5),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.0),
    minNrOfHitsForRebuild = cms.int32(2),
    propagatorAlong = cms.string('PropagatorWithMaterial'),
    propagatorOpposite = cms.string('PropagatorWithMaterialOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    seedAs5DHit = cms.bool(False),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('muonSeededTrajectoryFilterForInOut')
    ),
    updator = cms.string('KFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.muonSeededTrajectoryBuilderForOutIn = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    TTRHBuilder = cms.string('WithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('muonSeededMeasurementEstimatorForOutIn'),
    foundHitBonus = cms.double(1000.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('muonSeededTrajectoryFilterForOutIn')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(1.0),
    maxCand = cms.int32(3),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.0),
    minNrOfHitsForRebuild = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterial'),
    propagatorOpposite = cms.string('PropagatorWithMaterialOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    seedAs5DHit = cms.bool(False),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('muonSeededTrajectoryFilterForOutIn')
    ),
    updator = cms.string('KFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.muonSeededTrajectoryFilterForInOut = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(10),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    highEtaSwitch = cms.double(5.0),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(999),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutNone')
    ),
    minHitsAtHighEta = cms.int32(5),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.9),
    minimumNumberOfHits = cms.int32(3),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.muonSeededTrajectoryFilterForOutIn = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(10),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    highEtaSwitch = cms.double(5.0),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(999),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutNone')
    ),
    minHitsAtHighEta = cms.int32(5),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.9),
    minimumNumberOfHits = cms.int32(5),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.noPrediscriminants = cms.PSet(
    BooleanOperator = cms.string('and')
)

process.options = cms.untracked.PSet(
    FailPath = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    SkipEvent = cms.untracked.vstring(),
    accelerators = cms.untracked.vstring('*'),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    deleteNonConsumedUnscheduledModules = cms.untracked.bool(True),
    dumpOptions = cms.untracked.bool(False),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(0)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    holdsReferencesToDeleteEarly = cms.untracked.VPSet(),
    makeTriggerResults = cms.obsolete.untracked.bool,
    modulesToIgnoreForDeleteEarly = cms.untracked.vstring(),
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True)
)

process.pixelLessStepTrajectoryBuilder = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    TTRHBuilder = cms.string('WithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(False),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('pixelLessStepChi2Est'),
    foundHitBonus = cms.double(10.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('CkfBaseTrajectoryFilter_block')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(2),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.0),
    minNrOfHitsForRebuild = cms.int32(4),
    propagatorAlong = cms.string('PropagatorWithMaterial'),
    propagatorOpposite = cms.string('PropagatorWithMaterialOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    seedAs5DHit = cms.bool(False),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('pixelLessStepTrajectoryFilter')
    ),
    updator = cms.string('KFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.pixelLessStepTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    highEtaSwitch = cms.double(5.0),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(0),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutNone')
    ),
    minHitsAtHighEta = cms.int32(5),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.1),
    minimumNumberOfHits = cms.int32(4),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(1),
    strictSeedExtension = cms.bool(False)
)

process.pixelPairStepTrajectoryBuilder = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    TTRHBuilder = cms.string('WithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('pixelPairStepChi2Est'),
    foundHitBonus = cms.double(10.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('pixelPairStepTrajectoryFilterInOut')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(3),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.0),
    minNrOfHitsForRebuild = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterial'),
    propagatorOpposite = cms.string('PropagatorWithMaterialOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    seedAs5DHit = cms.bool(False),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('pixelPairStepTrajectoryFilter')
    ),
    updator = cms.string('KFUpdator'),
    useSameTrajFilter = cms.bool(False)
)

process.pixelPairStepTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('CompositeTrajectoryFilter'),
    filters = cms.VPSet(cms.PSet(
        refToPSet_ = cms.string('pixelPairStepTrajectoryFilterBase')
    ))
)

process.pixelPairStepTrajectoryFilterBase = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    highEtaSwitch = cms.double(5.0),
    maxCCCLostHits = cms.int32(0),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(999),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutLoose')
    ),
    minHitsAtHighEta = cms.int32(5),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.1),
    minimumNumberOfHits = cms.int32(4),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.pixelPairStepTrajectoryFilterInOut = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    highEtaSwitch = cms.double(5.0),
    maxCCCLostHits = cms.int32(0),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(999),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutLoose')
    ),
    minHitsAtHighEta = cms.int32(5),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.1),
    minimumNumberOfHits = cms.int32(4),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(1),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.pixelPairStepTrajectoryFilterShape = cms.PSet(
    ComponentType = cms.string('StripSubClusterShapeTrajectoryFilter'),
    layerMask = cms.PSet(
        TEC = cms.bool(False),
        TIB = cms.vuint32(1, 2),
        TID = cms.vuint32(1, 2),
        TOB = cms.bool(False)
    ),
    maxNSat = cms.uint32(3),
    maxTrimmedSizeDiffNeg = cms.double(1.0),
    maxTrimmedSizeDiffPos = cms.double(0.7),
    seedCutMIPs = cms.double(0.35),
    seedCutSN = cms.double(7.0),
    subclusterCutMIPs = cms.double(0.45),
    subclusterCutSN = cms.double(12.0),
    subclusterWindow = cms.double(0.7),
    trimMaxADC = cms.double(30.0),
    trimMaxFracNeigh = cms.double(0.25),
    trimMaxFracTotal = cms.double(0.15)
)

process.pset = cms.PSet(
    IDname = cms.string('ByTightCombinedIsolationDBSumPtCorr3HitsdR03'),
    maximumAbsoluteValues = cms.vdouble(0.8, 1000000000.0),
    maximumRelativeValues = cms.vdouble(-1.0, 0.1),
    referenceRawIDNames = cms.vstring(
        'ByRawCombinedIsolationDBSumPtCorr3HitsdR03',
        'PhotonPtSumOutsideSignalConedR03'
    )
)

process.qualityCutDictionary = cms.PSet(
    DetachedQuadStep = cms.vdouble(-0.66, -0.15, 0.46),
    DetachedTripletStep = cms.vdouble(-0.42, 0.16, 0.78),
    HighPtTripletStep = cms.vdouble(0.47, 0.55, 0.62),
    InitialStep = cms.vdouble(-0.35, 0.1, 0.28),
    JetCoreRegionalStep = cms.vdouble(-0.62, -0.49, 0.02),
    LowPtQuadStep = cms.vdouble(-0.37, 0.08, 0.28),
    LowPtTripletStep = cms.vdouble(-0.26, 0.09, 0.33),
    MixedTripletStep = cms.vdouble(-0.86, -0.68, -0.43),
    PixelLessStep = cms.vdouble(-0.8, -0.69, -0.4),
    PixelPairStep = cms.vdouble(-0.31, -0.13, 0.13),
    TobTecStep = cms.vdouble(-0.76, -0.65, -0.55)
)

process.requireDecayMode = cms.PSet(
    BooleanOperator = cms.string('and'),
    decayMode = cms.PSet(
        Producer = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs"),
        cut = cms.double(0.5)
    )
)

process.requireLeadPion = cms.PSet(
    BooleanOperator = cms.string('and'),
    leadPion = cms.PSet(
        Producer = cms.InputTag("pfRecoTauDiscriminationByLeadingTrackFinding"),
        cut = cms.double(0.5)
    )
)

process.requireLeadTrack = cms.PSet(
    BooleanOperator = cms.string('and'),
    leadTrack = cms.PSet(
        Producer = cms.InputTag("pfRecoTauDiscriminationByLeadingTrackFinding"),
        cut = cms.double(0.5)
    )
)

process.tautagInfoModifer = cms.PSet(
    name = cms.string('TTIworkaround'),
    pfTauTagInfoSrc = cms.InputTag("pfRecoTauTagInfoProducer"),
    plugin = cms.string('RecoTauTagInfoWorkaroundModifer')
)

process.tobTecStepInOutTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    highEtaSwitch = cms.double(5.0),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(0),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutNone')
    ),
    minHitsAtHighEta = cms.int32(5),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.1),
    minimumNumberOfHits = cms.int32(4),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(1),
    strictSeedExtension = cms.bool(False)
)

process.tobTecStepTrajectoryBuilder = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    TTRHBuilder = cms.string('WithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(False),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('tobTecStepChi2Est'),
    foundHitBonus = cms.double(10.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('tobTecStepInOutTrajectoryFilter')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(2),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.0),
    minNrOfHitsForRebuild = cms.int32(4),
    propagatorAlong = cms.string('PropagatorWithMaterial'),
    propagatorOpposite = cms.string('PropagatorWithMaterialOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    seedAs5DHit = cms.bool(False),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('tobTecStepTrajectoryFilter')
    ),
    updator = cms.string('KFUpdator'),
    useSameTrajFilter = cms.bool(False)
)

process.tobTecStepTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    highEtaSwitch = cms.double(5.0),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(0),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutNone')
    ),
    minHitsAtHighEta = cms.int32(5),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.1),
    minimumNumberOfHits = cms.int32(5),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(1),
    strictSeedExtension = cms.bool(False)
)

process.ClassifierMerger = cms.EDProducer("ClassifierMerger",
    inputClassifiers = cms.vstring(),
    mightGet = cms.optional.untracked.vstring
)


process.DuplicateListMerger = cms.EDProducer("DuplicateListMerger",
    candidateComponents = cms.InputTag(""),
    candidateSource = cms.InputTag(""),
    copyExtras = cms.untracked.bool(True),
    copyTrajectories = cms.untracked.bool(False),
    diffHitsCut = cms.int32(5),
    mergedMVAVals = cms.InputTag(""),
    mergedSource = cms.InputTag(""),
    mightGet = cms.optional.untracked.vstring,
    originalMVAVals = cms.InputTag(""),
    originalSource = cms.InputTag(""),
    trackAlgoPriorityOrder = cms.string('trackAlgoPriorityOrder')
)


process.DuplicateTrackMerger = cms.EDProducer("DuplicateTrackMerger",
    GBRForestFileName = cms.string(''),
    chi2EstimatorName = cms.string('DuplicateTrackMergerChi2Est'),
    forestLabel = cms.string('MVADuplicate'),
    maxDCA = cms.double(30),
    maxDLambda = cms.double(0.3),
    maxDPhi = cms.double(0.3),
    maxDQoP = cms.double(0.25),
    maxDdsz = cms.double(10),
    maxDdxy = cms.double(10),
    mightGet = cms.optional.untracked.vstring,
    minBDTG = cms.double(-0.1),
    minDeltaR3d = cms.double(-4),
    minP = cms.double(0.4),
    minpT = cms.double(0.2),
    overlapCheckMaxHits = cms.uint32(4),
    overlapCheckMaxMissingLayers = cms.uint32(1),
    overlapCheckMinCosT = cms.double(0.99),
    propagatorName = cms.string('PropagatorWithMaterial'),
    source = cms.InputTag(""),
    ttrhBuilderName = cms.string('WithAngleAndTemplate'),
    useInnermostState = cms.bool(True)
)


process.MeasurementTrackerEvent = cms.EDProducer("MeasurementTrackerEventProducer",
    Phase2TrackerCluster1DProducer = cms.string(''),
    badPixelFEDChannelCollectionLabels = cms.VInputTag("siPixelDigis"),
    inactivePixelDetectorLabels = cms.VInputTag(),
    inactiveStripDetectorLabels = cms.VInputTag("siStripDigis"),
    measurementTracker = cms.string(''),
    mightGet = cms.optional.untracked.vstring,
    pixelCablingMapLabel = cms.string(''),
    pixelClusterProducer = cms.string('siPixelClusters'),
    skipClusters = cms.InputTag(""),
    stripClusterProducer = cms.string('siStripClusters'),
    switchOffPixelsIfEmpty = cms.bool(True),
    vectorHits = cms.InputTag(""),
    vectorHitsRej = cms.InputTag("")
)


process.MeasurementTrackerEventPreSplitting = cms.EDProducer("MeasurementTrackerEventProducer",
    Phase2TrackerCluster1DProducer = cms.string(''),
    badPixelFEDChannelCollectionLabels = cms.VInputTag("siPixelDigis"),
    inactivePixelDetectorLabels = cms.VInputTag("siPixelDigis"),
    inactiveStripDetectorLabels = cms.VInputTag("siStripDigis"),
    measurementTracker = cms.string(''),
    mightGet = cms.optional.untracked.vstring,
    pixelCablingMapLabel = cms.string(''),
    pixelClusterProducer = cms.string('siPixelClustersPreSplitting'),
    skipClusters = cms.InputTag(""),
    stripClusterProducer = cms.string('siStripClusters'),
    switchOffPixelsIfEmpty = cms.bool(True),
    vectorHits = cms.InputTag(""),
    vectorHitsRej = cms.InputTag("")
)


process.MixedLayerPairs = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(
        HitProducer = cms.string('siPixelRecHits'),
        TTRHBuilder = cms.string('WithTrackAngle')
    ),
    FPix = cms.PSet(
        HitProducer = cms.string('siPixelRecHits'),
        TTRHBuilder = cms.string('WithTrackAngle')
    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        maxRing = cms.int32(1),
        minRing = cms.int32(1),
        useRingSlector = cms.bool(True)
    ),
    TIB = cms.PSet(

    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(

    ),
    layerList = cms.vstring(
        'BPix1+BPix2',
        'BPix1+BPix3',
        'BPix2+BPix3',
        'BPix1+FPix1_pos',
        'BPix1+FPix1_neg',
        'BPix1+FPix2_pos',
        'BPix1+FPix2_neg',
        'BPix2+FPix1_pos',
        'BPix2+FPix1_neg',
        'BPix2+FPix2_pos',
        'BPix2+FPix2_neg',
        'FPix1_pos+FPix2_pos',
        'FPix1_neg+FPix2_neg',
        'FPix2_pos+TEC1_pos',
        'FPix2_pos+TEC2_pos',
        'TEC1_pos+TEC2_pos',
        'TEC2_pos+TEC3_pos',
        'FPix2_neg+TEC1_neg',
        'FPix2_neg+TEC2_neg',
        'TEC1_neg+TEC2_neg',
        'TEC2_neg+TEC3_neg'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.MixedLayerTriplets = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(
        HitProducer = cms.string('siPixelRecHits'),
        TTRHBuilder = cms.string('WithTrackAngle')
    ),
    FPix = cms.PSet(
        HitProducer = cms.string('siPixelRecHits'),
        TTRHBuilder = cms.string('WithTrackAngle')
    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit")
    ),
    TIB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit")
    ),
    TID = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit")
    ),
    TOB = cms.PSet(

    ),
    layerList = cms.vstring(
        'BPix1+BPix2+BPix3',
        'BPix1+BPix2+FPix1_pos',
        'BPix1+BPix2+FPix1_neg',
        'BPix1+FPix1_pos+FPix2_pos',
        'BPix1+FPix1_neg+FPix2_neg',
        'BPix1+BPix2+TIB1',
        'BPix1+BPix3+TIB1',
        'BPix2+BPix3+TIB1',
        'BPix1+FPix1_pos+TID1_pos',
        'BPix1+FPix1_neg+TID1_neg',
        'BPix1+FPix1_pos+TID2_pos',
        'BPix1+FPix1_neg+TID2_neg',
        'FPix1_pos+FPix2_pos+TEC1_pos',
        'FPix1_neg+FPix2_neg+TEC1_neg',
        'FPix1_pos+FPix2_pos+TEC2_pos',
        'FPix1_neg+FPix2_neg+TEC2_neg'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.PFTauPrimaryVertexProducer = cms.EDProducer("PFTauPrimaryVertexProducer",
    Algorithm = cms.int32(0),
    ElectronTag = cms.InputTag("MyElectrons"),
    MuonTag = cms.InputTag("MyMuons"),
    PFTauTag = cms.InputTag("hpsPFTauProducer"),
    PVTag = cms.InputTag("offlinePrimaryVertices"),
    RemoveElectronTracks = cms.bool(False),
    RemoveMuonTracks = cms.bool(False),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    cut = cms.string('pt > 18.0 & abs(eta)<2.3'),
    discriminators = cms.VPSet(cms.PSet(
        discriminator = cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding"),
        selectionCut = cms.double(0.5)
    )),
    mightGet = cms.optional.untracked.vstring,
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.2),
            maxDeltaZToLeadTrack = cms.double(-1.0),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.03),
            minGammaEt = cms.double(1.5),
            minTrackHits = cms.uint32(8),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(1.0),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        leadingTrkOrPFCandOption = cms.string('leadPFCand'),
        primaryVertexSrc = cms.InputTag("offlinePrimaryVertices"),
        pvFindingAlgo = cms.string('closestInDeltaZ'),
        recoverLeadingTrk = cms.bool(False),
        signalQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.4),
            maxDeltaZToLeadTrack = cms.double(-1.0),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(1.0),
            minNeutralHadronEt = cms.double(30.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        vertexTrackFiltering = cms.bool(False),
        vxAssocQualityCuts = cms.PSet(
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(1.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        )
    ),
    useBeamSpot = cms.bool(True),
    useSelectedTaus = cms.bool(False)
)


process.PFTauSecondaryVertexProducer = cms.EDProducer("PFTauSecondaryVertexProducer",
    PFTauTag = cms.InputTag("hpsPFTauProducer"),
    mightGet = cms.optional.untracked.vstring
)


process.PFTauTransverseImpactParameters = cms.EDProducer("PFTauTransverseImpactParameters",
    PFTauPVATag = cms.InputTag("PFTauPrimaryVertexProducer"),
    PFTauSVATag = cms.InputTag("PFTauSecondaryVertexProducer"),
    PFTauTag = cms.InputTag("hpsPFTauProducer"),
    mightGet = cms.optional.untracked.vstring,
    useFullCalculation = cms.bool(False)
)


process.PixelLayerPairs = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(
        HitProducer = cms.string('siPixelRecHits'),
        TTRHBuilder = cms.string('WithTrackAngle')
    ),
    FPix = cms.PSet(
        HitProducer = cms.string('siPixelRecHits'),
        TTRHBuilder = cms.string('WithTrackAngle')
    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(

    ),
    TIB = cms.PSet(

    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(

    ),
    layerList = cms.vstring(
        'BPix1+BPix2',
        'BPix1+BPix3',
        'BPix2+BPix3',
        'BPix1+FPix1_pos',
        'BPix1+FPix1_neg',
        'BPix1+FPix2_pos',
        'BPix1+FPix2_neg',
        'BPix2+FPix1_pos',
        'BPix2+FPix1_neg',
        'BPix2+FPix2_pos',
        'BPix2+FPix2_neg',
        'FPix1_pos+FPix2_pos',
        'FPix1_neg+FPix2_neg'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.PixelLayerTriplets = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(
        HitProducer = cms.string('siPixelRecHits'),
        TTRHBuilder = cms.string('WithTrackAngle')
    ),
    FPix = cms.PSet(
        HitProducer = cms.string('siPixelRecHits'),
        TTRHBuilder = cms.string('WithTrackAngle')
    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(

    ),
    TIB = cms.PSet(

    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(

    ),
    layerList = cms.vstring(
        'BPix1+BPix2+BPix3',
        'BPix2+BPix3+BPix4',
        'BPix1+BPix3+BPix4',
        'BPix1+BPix2+BPix4',
        'BPix2+BPix3+FPix1_pos',
        'BPix2+BPix3+FPix1_neg',
        'BPix1+BPix2+FPix1_pos',
        'BPix1+BPix2+FPix1_neg',
        'BPix2+FPix1_pos+FPix2_pos',
        'BPix2+FPix1_neg+FPix2_neg',
        'BPix1+FPix1_pos+FPix2_pos',
        'BPix1+FPix1_neg+FPix2_neg',
        'FPix1_pos+FPix2_pos+FPix3_pos',
        'FPix1_neg+FPix2_neg+FPix3_neg'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.RecoTauCleaner = cms.EDProducer("RecoTauCleaner",
    cleaners = cms.VPSet(
        cms.PSet(
            name = cms.string('Charge'),
            nprongs = cms.vuint32(1, 3),
            passForCharge = cms.int32(1),
            plugin = cms.string('RecoTauChargeCleanerPlugin'),
            selectionFailValue = cms.double(0),
            tolerance = cms.double(0)
        ),
        cms.PSet(
            name = cms.string('HPS_Select'),
            plugin = cms.string('RecoTauDiscriminantCleanerPlugin'),
            src = cms.InputTag("hpsSelectionDiscriminator"),
            tolerance = cms.double(0)
        ),
        cms.PSet(
            minTrackPt = cms.double(5.0),
            name = cms.string('killSoftTwoProngTaus'),
            plugin = cms.string('RecoTauSoftTwoProngTausCleanerPlugin'),
            tolerance = cms.double(0)
        ),
        cms.PSet(
            name = cms.string('ChargedHadronMultiplicity'),
            plugin = cms.string('RecoTauChargedHadronMultiplicityCleanerPlugin'),
            tolerance = cms.double(0)
        ),
        cms.PSet(
            name = cms.string('Pt'),
            plugin = cms.string('RecoTauStringCleanerPlugin'),
            selection = cms.string('leadCand().isNonnull()'),
            selectionFailValue = cms.double(1000.0),
            selectionPassFunction = cms.string('-pt()'),
            tolerance = cms.double(0.01)
        ),
        cms.PSet(
            name = cms.string('StripMultiplicity'),
            plugin = cms.string('RecoTauStringCleanerPlugin'),
            selection = cms.string('leadCand().isNonnull()'),
            selectionFailValue = cms.double(1000.0),
            selectionPassFunction = cms.string('-signalPiZeroCandidates().size()'),
            tolerance = cms.double(0)
        ),
        cms.PSet(
            name = cms.string('CombinedIsolation'),
            plugin = cms.string('RecoTauStringCleanerPlugin'),
            selection = cms.string('leadCand().isNonnull()'),
            selectionFailValue = cms.double(1000.0),
            selectionPassFunction = cms.string('isolationPFChargedHadrCandsPtSum() + isolationPFGammaCandsEtSum()'),
            tolerance = cms.double(0)
        )
    ),
    outputSelection = cms.string(''),
    src = cms.InputTag("combinatoricRecoTaus"),
    verbosity = cms.int32(0)
)


process.RecoTauJetRegionProducer = cms.EDProducer("RecoTauJetRegionProducer",
    deltaR = cms.double(0.8),
    maxJetAbsEta = cms.double(2.7),
    minJetPt = cms.double(14.0),
    pfCandAssocMapSrc = cms.InputTag(""),
    pfCandSrc = cms.InputTag("particleFlow"),
    src = cms.InputTag("ak4PFJets"),
    verbosity = cms.int32(0)
)


process.RecoTauPiZeroUnembedder = cms.EDProducer("RecoTauPiZeroUnembedder",
    src = cms.InputTag("hpsPFTauProducerSansRefs")
)


process.SiStripApprox2Clusters = cms.EDProducer("SiStripApprox2Clusters",
    inputApproxClusters = cms.InputTag("SiStripClusters2ApproxClusters"),
    mightGet = cms.optional.untracked.vstring
)


process.TrackCollectionMerger = cms.EDProducer("TrackCollectionMerger",
    allowFirstHitShare = cms.bool(True),
    copyExtras = cms.untracked.bool(True),
    copyTrajectories = cms.untracked.bool(False),
    enableMerging = cms.bool(True),
    foundHitBonus = cms.double(10),
    inputClassifiers = cms.vstring(),
    lostHitPenalty = cms.double(5),
    mightGet = cms.optional.untracked.vstring,
    minQuality = cms.string('loose'),
    minShareHits = cms.uint32(2),
    shareFrac = cms.double(0.19),
    trackAlgoPriorityOrder = cms.string('trackAlgoPriorityOrder'),
    trackProducers = cms.VInputTag()
)


process.TrackCutClassifier = cms.EDProducer("TrackCutClassifier",
    beamspot = cms.InputTag("offlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mightGet = cms.optional.untracked.vstring,
    mva = cms.PSet(
        dr_par = cms.PSet(
            d0err = cms.vdouble(0.003, 0.003, 0.003),
            d0err_par = cms.vdouble(0.001, 0.001, 0.001),
            drWPVerr_par = cms.vdouble(3.4028234663852886e+38, 3.4028234663852886e+38, 3.4028234663852886e+38),
            dr_exp = cms.vint32(2147483647, 2147483647, 2147483647),
            dr_par1 = cms.vdouble(3.4028234663852886e+38, 3.4028234663852886e+38, 3.4028234663852886e+38),
            dr_par2 = cms.vdouble(3.4028234663852886e+38, 3.4028234663852886e+38, 3.4028234663852886e+38)
        ),
        dz_par = cms.PSet(
            dzWPVerr_par = cms.vdouble(3.4028234663852886e+38, 3.4028234663852886e+38, 3.4028234663852886e+38),
            dz_exp = cms.vint32(2147483647, 2147483647, 2147483647),
            dz_par1 = cms.vdouble(3.4028234663852886e+38, 3.4028234663852886e+38, 3.4028234663852886e+38),
            dz_par2 = cms.vdouble(3.4028234663852886e+38, 3.4028234663852886e+38, 3.4028234663852886e+38)
        ),
        isHLT = cms.bool(False),
        maxChi2 = cms.vdouble(9999, 25, 16),
        maxChi2n = cms.vdouble(9999, 1, 0.4),
        maxDr = cms.vdouble(3.4028234663852886e+38, 3.4028234663852886e+38, 3.4028234663852886e+38),
        maxDz = cms.vdouble(3.4028234663852886e+38, 3.4028234663852886e+38, 3.4028234663852886e+38),
        maxDzWrtBS = cms.vdouble(3.4028234663852886e+38, 24, 15),
        maxLostLayers = cms.vint32(99, 3, 3),
        maxRelPtErr = cms.vdouble(3.4028234663852886e+38, 3.4028234663852886e+38, 3.4028234663852886e+38),
        min3DLayers = cms.vint32(1, 2, 3),
        minHits = cms.vint32(0, 0, 1),
        minHits4pass = cms.vint32(2147483647, 2147483647, 2147483647),
        minLayers = cms.vint32(3, 4, 5),
        minNVtxTrk = cms.int32(2),
        minNdof = cms.vdouble(-1, -1, -1),
        minPixelHits = cms.vint32(0, 0, 1)
    ),
    qualityCuts = cms.vdouble(-0.7, 0.1, 0.7),
    src = cms.InputTag(""),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.TrackMVAClassifierDetached = cms.EDProducer("TrackMVAClassifierDetached",
    beamspot = cms.InputTag("offlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mightGet = cms.optional.untracked.vstring,
    mva = cms.PSet(
        GBRForestFileName = cms.string(''),
        GBRForestLabel = cms.string('')
    ),
    qualityCuts = cms.vdouble(-0.7, 0.1, 0.7),
    src = cms.InputTag(""),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.TrackMVAClassifierPrompt = cms.EDProducer("TrackMVAClassifierPrompt",
    beamspot = cms.InputTag("offlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mightGet = cms.optional.untracked.vstring,
    mva = cms.PSet(
        GBRForestFileName = cms.string(''),
        GBRForestLabel = cms.string('')
    ),
    qualityCuts = cms.vdouble(-0.7, 0.1, 0.7),
    src = cms.InputTag(""),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.TrackProducer = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('undefAlgorithm'),
    Fitter = cms.string('KFFittingSmootherWithOutliersRejectionAndRK'),
    GeometricInnerState = cms.bool(False),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    Propagator = cms.string('RungeKuttaTrackerPropagator'),
    SimpleMagneticField = cms.string(''),
    TTRHBuilder = cms.string('WithAngleAndTemplate'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("ckfTrackCandidates"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(False)
)


process.TrackRefitter = cms.EDProducer("TrackRefitter",
    AlgorithmName = cms.string('undefAlgorithm'),
    Fitter = cms.string('KFFittingSmootherWithOutliersRejectionAndRK'),
    GeometricInnerState = cms.bool(False),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    Propagator = cms.string('RungeKuttaTrackerPropagator'),
    TTRHBuilder = cms.string('WithAngleAndTemplate'),
    TrajectoryInEvent = cms.bool(True),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    constraint = cms.string(''),
    src = cms.InputTag("generalTracks"),
    srcConstr = cms.InputTag(""),
    useHitsSplitting = cms.bool(False)
)


process.TrackRefitterBHM = cms.EDProducer("TrackRefitter",
    AlgorithmName = cms.string('undefAlgorithm'),
    Fitter = cms.string('KFFittingSmootherBH'),
    GeometricInnerState = cms.bool(True),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    Propagator = cms.string('BeamHaloPropagatorAlong'),
    TTRHBuilder = cms.string('WithAngleAndTemplate'),
    TrajectoryInEvent = cms.bool(True),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    constraint = cms.string(''),
    src = cms.InputTag("ctfWithMaterialTracksBeamHaloMuon"),
    srcConstr = cms.InputTag(""),
    useHitsSplitting = cms.bool(False)
)


process.TrackRefitterP5 = cms.EDProducer("TrackRefitter",
    AlgorithmName = cms.string('ctf'),
    Fitter = cms.string('FittingSmootherRKP5'),
    GeometricInnerState = cms.bool(True),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    Propagator = cms.string('RungeKuttaTrackerPropagator'),
    TTRHBuilder = cms.string('WithAngleAndTemplate'),
    TrajectoryInEvent = cms.bool(True),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    constraint = cms.string(''),
    src = cms.InputTag("ctfWithMaterialTracksP5"),
    srcConstr = cms.InputTag(""),
    useHitsSplitting = cms.bool(False)
)


process.ak4CaloJetsForTrk = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(1),
    GhostArea = cms.double(0.01),
    Ghost_EtaMax = cms.double(5.0),
    Rho_EtaMax = cms.double(4.4),
    applyWeight = cms.bool(False),
    doAreaDiskApprox = cms.bool(False),
    doAreaFastjet = cms.bool(False),
    doPUOffsetCorr = cms.bool(False),
    doPVCorrection = cms.bool(True),
    doRhoFastjet = cms.bool(False),
    inputEMin = cms.double(0.0),
    inputEtMin = cms.double(0.3),
    jetAlgorithm = cms.string('AntiKt'),
    jetPtMin = cms.double(10.0),
    jetType = cms.string('CaloJet'),
    maxBadEcalCells = cms.uint32(9999999),
    maxBadHcalCells = cms.uint32(9999999),
    maxProblematicEcalCells = cms.uint32(9999999),
    maxProblematicHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    maxRecoveredHcalCells = cms.uint32(9999999),
    minSeed = cms.uint32(14327),
    nSigmaPU = cms.double(1.0),
    puPtMin = cms.double(10),
    rParam = cms.double(0.4),
    radiusPU = cms.double(0.5),
    src = cms.InputTag("caloTowerForTrk"),
    srcPVs = cms.InputTag("firstStepPrimaryVerticesUnsorted"),
    useDeterministicSeed = cms.bool(True),
    voronoiRfact = cms.double(-0.9)
)


process.ak4CaloJetsForTrkPreSplitting = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(1),
    GhostArea = cms.double(0.01),
    Ghost_EtaMax = cms.double(5.0),
    Rho_EtaMax = cms.double(4.4),
    applyWeight = cms.bool(False),
    doAreaDiskApprox = cms.bool(False),
    doAreaFastjet = cms.bool(False),
    doPUOffsetCorr = cms.bool(False),
    doPVCorrection = cms.bool(True),
    doRhoFastjet = cms.bool(False),
    inputEMin = cms.double(0.0),
    inputEtMin = cms.double(0.3),
    jetAlgorithm = cms.string('AntiKt'),
    jetPtMin = cms.double(10.0),
    jetType = cms.string('CaloJet'),
    maxBadEcalCells = cms.uint32(9999999),
    maxBadHcalCells = cms.uint32(9999999),
    maxProblematicEcalCells = cms.uint32(9999999),
    maxProblematicHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    maxRecoveredHcalCells = cms.uint32(9999999),
    minSeed = cms.uint32(14327),
    nSigmaPU = cms.double(1.0),
    puPtMin = cms.double(10),
    rParam = cms.double(0.4),
    radiusPU = cms.double(0.5),
    src = cms.InputTag("caloTowerForTrkPreSplitting"),
    srcPVs = cms.InputTag("firstStepPrimaryVerticesPreSplitting"),
    useDeterministicSeed = cms.bool(True),
    voronoiRfact = cms.double(-0.9)
)


process.ak4JetTracksAssociatorAtVertexPF = cms.EDProducer("JetTracksAssociatorAtVertex",
    coneSize = cms.double(0.4),
    jets = cms.InputTag("ak4PFJetsCHS"),
    pvSrc = cms.InputTag("offlinePrimaryVertices"),
    tracks = cms.InputTag("generalTracks"),
    useAssigned = cms.bool(False)
)


process.ak4PFJetTracksAssociatorAtVertex = cms.EDProducer("JetTracksAssociatorAtVertex",
    coneSize = cms.double(0.4),
    jets = cms.InputTag("ak4PFJets"),
    pvSrc = cms.InputTag("offlinePrimaryVertices"),
    tracks = cms.InputTag("generalTracks"),
    useAssigned = cms.bool(False)
)


process.ak4PFJetsLegacyHPSPiZeros = cms.EDProducer("RecoTauPiZeroProducer",
    builders = cms.VPSet(cms.PSet(
        applyElecTrackQcuts = cms.bool(False),
        makeCombinatoricStrips = cms.bool(False),
        maxStripBuildIterations = cms.int32(-1),
        minGammaEtStripAdd = cms.double(1.0),
        minGammaEtStripSeed = cms.double(1.0),
        minStripEt = cms.double(1.0),
        name = cms.string('s'),
        plugin = cms.string('RecoTauPiZeroStripPlugin3'),
        qualityCuts = cms.PSet(
            isolationQualityCuts = cms.PSet(
                maxDeltaZ = cms.double(0.2),
                maxDeltaZToLeadTrack = cms.double(-1.0),
                maxTrackChi2 = cms.double(100.0),
                maxTransverseImpactParameter = cms.double(0.03),
                minGammaEt = cms.double(1.5),
                minTrackHits = cms.uint32(8),
                minTrackPixelHits = cms.uint32(0),
                minTrackPt = cms.double(1.0),
                minTrackVertexWeight = cms.double(-1.0)
            ),
            leadingTrkOrPFCandOption = cms.string('leadPFCand'),
            primaryVertexSrc = cms.InputTag("offlinePrimaryVertices"),
            pvFindingAlgo = cms.string('closestInDeltaZ'),
            recoverLeadingTrk = cms.bool(False),
            signalQualityCuts = cms.PSet(
                maxDeltaZ = cms.double(0.4),
                maxDeltaZToLeadTrack = cms.double(-1.0),
                maxTrackChi2 = cms.double(100.0),
                maxTransverseImpactParameter = cms.double(0.1),
                minGammaEt = cms.double(1.0),
                minNeutralHadronEt = cms.double(30.0),
                minTrackHits = cms.uint32(3),
                minTrackPixelHits = cms.uint32(0),
                minTrackPt = cms.double(0.5),
                minTrackVertexWeight = cms.double(-1.0)
            ),
            vertexTrackFiltering = cms.bool(False),
            vxAssocQualityCuts = cms.PSet(
                maxTrackChi2 = cms.double(100.0),
                maxTransverseImpactParameter = cms.double(0.1),
                minGammaEt = cms.double(1.0),
                minTrackHits = cms.uint32(3),
                minTrackPixelHits = cms.uint32(0),
                minTrackPt = cms.double(0.5),
                minTrackVertexWeight = cms.double(-1.0)
            )
        ),
        stripCandidatesParticleIds = cms.vint32(2, 4),
        stripEtaAssociationDistance = cms.double(0.05),
        stripEtaAssociationDistanceFunc = cms.PSet(
            function = cms.string('TMath::Min(0.15, TMath::Max(0.05, [0]*TMath::Power(pT, -[1])))'),
            par0 = cms.double(0.197077),
            par1 = cms.double(0.658701)
        ),
        stripPhiAssociationDistance = cms.double(0.2),
        stripPhiAssociationDistanceFunc = cms.PSet(
            function = cms.string('TMath::Min(0.3, TMath::Max(0.05, [0]*TMath::Power(pT, -[1])))'),
            par0 = cms.double(0.352476),
            par1 = cms.double(0.707716)
        ),
        updateStripAfterEachDaughter = cms.bool(False),
        verbosity = cms.int32(0)
    )),
    jetSrc = cms.InputTag("ak4PFJets"),
    massHypothesis = cms.double(0.136),
    maxJetAbsEta = cms.double(2.7),
    mightGet = cms.optional.untracked.vstring,
    minJetPt = cms.double(14.0),
    outputSelection = cms.string('pt > 0'),
    ranking = cms.VPSet(cms.PSet(
        name = cms.string('InStrip'),
        plugin = cms.string('RecoTauPiZeroStringQuality'),
        selection = cms.string('algoIs("kStrips")'),
        selectionFailValue = cms.double(1000),
        selectionPassFunction = cms.string('abs(mass() - 0.13579)')
    )),
    verbosity = cms.int32(0)
)


process.ak4PFJetsRecoTauChargedHadrons = cms.EDProducer("PFRecoTauChargedHadronProducer",
    builders = cms.VPSet(
        cms.PSet(
            chargedHadronCandidatesParticleIds = cms.vint32(1, 2, 3),
            dRmergeNeutralHadronWrtChargedHadron = cms.double(0.005),
            dRmergeNeutralHadronWrtElectron = cms.double(0.05),
            dRmergeNeutralHadronWrtNeutralHadron = cms.double(0.01),
            dRmergeNeutralHadronWrtOther = cms.double(0.005),
            dRmergePhotonWrtChargedHadron = cms.double(0.005),
            dRmergePhotonWrtElectron = cms.double(0.005),
            dRmergePhotonWrtNeutralHadron = cms.double(0.01),
            dRmergePhotonWrtOther = cms.double(0.005),
            maxUnmatchedBlockElementsNeutralHadron = cms.int32(1),
            maxUnmatchedBlockElementsPhoton = cms.int32(1),
            minBlockElementMatchesNeutralHadron = cms.int32(2),
            minBlockElementMatchesPhoton = cms.int32(2),
            minMergeChargedHadronPt = cms.double(100.0),
            minMergeGammaEt = cms.double(1.0),
            minMergeNeutralHadronEt = cms.double(1.0),
            name = cms.string('chargedPFCandidates'),
            plugin = cms.string('PFRecoTauChargedHadronFromPFCandidatePlugin'),
            qualityCuts = cms.PSet(
                isolationQualityCuts = cms.PSet(
                    maxDeltaZ = cms.double(0.2),
                    maxDeltaZToLeadTrack = cms.double(-1.0),
                    maxTrackChi2 = cms.double(100.0),
                    maxTransverseImpactParameter = cms.double(0.03),
                    minGammaEt = cms.double(1.5),
                    minTrackHits = cms.uint32(8),
                    minTrackPixelHits = cms.uint32(0),
                    minTrackPt = cms.double(1.0),
                    minTrackVertexWeight = cms.double(-1.0)
                ),
                leadingTrkOrPFCandOption = cms.string('leadPFCand'),
                primaryVertexSrc = cms.InputTag("offlinePrimaryVertices"),
                pvFindingAlgo = cms.string('closestInDeltaZ'),
                recoverLeadingTrk = cms.bool(False),
                signalQualityCuts = cms.PSet(
                    maxDeltaZ = cms.double(0.4),
                    maxDeltaZToLeadTrack = cms.double(-1.0),
                    maxTrackChi2 = cms.double(100.0),
                    maxTransverseImpactParameter = cms.double(0.1),
                    minGammaEt = cms.double(1.0),
                    minNeutralHadronEt = cms.double(30.0),
                    minTrackHits = cms.uint32(3),
                    minTrackPixelHits = cms.uint32(0),
                    minTrackPt = cms.double(0.5),
                    minTrackVertexWeight = cms.double(-1.0)
                ),
                vertexTrackFiltering = cms.bool(False),
                vxAssocQualityCuts = cms.PSet(
                    maxTrackChi2 = cms.double(100.0),
                    maxTransverseImpactParameter = cms.double(0.1),
                    minGammaEt = cms.double(1.0),
                    minTrackHits = cms.uint32(3),
                    minTrackPixelHits = cms.uint32(0),
                    minTrackPt = cms.double(0.5),
                    minTrackVertexWeight = cms.double(-1.0)
                )
            ),
            verbosity = cms.int32(0)
        ),
        cms.PSet(
            dRcone = cms.double(0.5),
            dRconeLimitedToJetArea = cms.bool(False),
            dRmergeNeutralHadron = cms.double(0.1),
            dRmergePhoton = cms.double(0.05),
            minMergeChargedHadronPt = cms.double(100.0),
            minMergeGammaEt = cms.double(1.0),
            minMergeNeutralHadronEt = cms.double(1.0),
            name = cms.string('tracks'),
            plugin = cms.string('PFRecoTauChargedHadronFromTrackPlugin'),
            qualityCuts = cms.PSet(
                isolationQualityCuts = cms.PSet(
                    maxDeltaZ = cms.double(0.2),
                    maxDeltaZToLeadTrack = cms.double(-1.0),
                    maxTrackChi2 = cms.double(100.0),
                    maxTransverseImpactParameter = cms.double(0.03),
                    minGammaEt = cms.double(1.5),
                    minTrackHits = cms.uint32(8),
                    minTrackPixelHits = cms.uint32(0),
                    minTrackPt = cms.double(1.0),
                    minTrackVertexWeight = cms.double(-1.0)
                ),
                leadingTrkOrPFCandOption = cms.string('leadPFCand'),
                primaryVertexSrc = cms.InputTag("offlinePrimaryVertices"),
                pvFindingAlgo = cms.string('closestInDeltaZ'),
                recoverLeadingTrk = cms.bool(False),
                signalQualityCuts = cms.PSet(
                    maxDeltaZ = cms.double(0.4),
                    maxDeltaZToLeadTrack = cms.double(-1.0),
                    maxTrackChi2 = cms.double(100.0),
                    maxTransverseImpactParameter = cms.double(0.1),
                    minGammaEt = cms.double(1.0),
                    minNeutralHadronEt = cms.double(30.0),
                    minTrackHits = cms.uint32(3),
                    minTrackPixelHits = cms.uint32(0),
                    minTrackPt = cms.double(0.5),
                    minTrackVertexWeight = cms.double(-1.0)
                ),
                vertexTrackFiltering = cms.bool(False),
                vxAssocQualityCuts = cms.PSet(
                    maxTrackChi2 = cms.double(100.0),
                    maxTransverseImpactParameter = cms.double(0.1),
                    minGammaEt = cms.double(1.0),
                    minTrackHits = cms.uint32(3),
                    minTrackPixelHits = cms.uint32(0),
                    minTrackPt = cms.double(0.5),
                    minTrackVertexWeight = cms.double(-1.0)
                )
            ),
            srcTracks = cms.InputTag("generalTracks"),
            verbosity = cms.int32(0)
        ),
        cms.PSet(
            chargedHadronCandidatesParticleIds = cms.vint32(5),
            dRmergeNeutralHadronWrtChargedHadron = cms.double(0.005),
            dRmergeNeutralHadronWrtElectron = cms.double(0.05),
            dRmergeNeutralHadronWrtNeutralHadron = cms.double(0.01),
            dRmergeNeutralHadronWrtOther = cms.double(0.005),
            dRmergePhotonWrtChargedHadron = cms.double(0.005),
            dRmergePhotonWrtElectron = cms.double(0.005),
            dRmergePhotonWrtNeutralHadron = cms.double(0.01),
            dRmergePhotonWrtOther = cms.double(0.005),
            maxUnmatchedBlockElementsNeutralHadron = cms.int32(1),
            maxUnmatchedBlockElementsPhoton = cms.int32(1),
            minBlockElementMatchesNeutralHadron = cms.int32(2),
            minBlockElementMatchesPhoton = cms.int32(2),
            minMergeChargedHadronPt = cms.double(0.0),
            minMergeGammaEt = cms.double(1.0),
            minMergeNeutralHadronEt = cms.double(1.0),
            name = cms.string('PFNeutralHadrons'),
            plugin = cms.string('PFRecoTauChargedHadronFromPFCandidatePlugin'),
            qualityCuts = cms.PSet(
                isolationQualityCuts = cms.PSet(
                    maxDeltaZ = cms.double(0.2),
                    maxDeltaZToLeadTrack = cms.double(-1.0),
                    maxTrackChi2 = cms.double(100.0),
                    maxTransverseImpactParameter = cms.double(0.03),
                    minGammaEt = cms.double(1.5),
                    minTrackHits = cms.uint32(8),
                    minTrackPixelHits = cms.uint32(0),
                    minTrackPt = cms.double(1.0),
                    minTrackVertexWeight = cms.double(-1.0)
                ),
                leadingTrkOrPFCandOption = cms.string('leadPFCand'),
                primaryVertexSrc = cms.InputTag("offlinePrimaryVertices"),
                pvFindingAlgo = cms.string('closestInDeltaZ'),
                recoverLeadingTrk = cms.bool(False),
                signalQualityCuts = cms.PSet(
                    maxDeltaZ = cms.double(0.4),
                    maxDeltaZToLeadTrack = cms.double(-1.0),
                    maxTrackChi2 = cms.double(100.0),
                    maxTransverseImpactParameter = cms.double(0.1),
                    minGammaEt = cms.double(1.0),
                    minNeutralHadronEt = cms.double(30.0),
                    minTrackHits = cms.uint32(3),
                    minTrackPixelHits = cms.uint32(0),
                    minTrackPt = cms.double(0.5),
                    minTrackVertexWeight = cms.double(-1.0)
                ),
                vertexTrackFiltering = cms.bool(False),
                vxAssocQualityCuts = cms.PSet(
                    maxTrackChi2 = cms.double(100.0),
                    maxTransverseImpactParameter = cms.double(0.1),
                    minGammaEt = cms.double(1.0),
                    minTrackHits = cms.uint32(3),
                    minTrackPixelHits = cms.uint32(0),
                    minTrackPt = cms.double(0.5),
                    minTrackVertexWeight = cms.double(-1.0)
                )
            ),
            verbosity = cms.int32(0)
        )
    ),
    jetSrc = cms.InputTag("ak4PFJets"),
    maxJetAbsEta = cms.double(2.7),
    mightGet = cms.optional.untracked.vstring,
    minJetPt = cms.double(14.0),
    outputSelection = cms.string('pt > 0.5'),
    ranking = cms.VPSet(
        cms.PSet(
            name = cms.string('ChargedPFCandidate'),
            plugin = cms.string('PFRecoTauChargedHadronStringQuality'),
            selection = cms.string("algoIs(\'kChargedPFCandidate\')"),
            selectionFailValue = cms.double(1000.0),
            selectionPassFunction = cms.string('-pt')
        ),
        cms.PSet(
            name = cms.string('ChargedPFCandidate'),
            plugin = cms.string('PFRecoTauChargedHadronStringQuality'),
            selection = cms.string("algoIs(\'kTrack\')"),
            selectionFailValue = cms.double(1000.0),
            selectionPassFunction = cms.string('-pt')
        ),
        cms.PSet(
            name = cms.string('ChargedPFCandidate'),
            plugin = cms.string('PFRecoTauChargedHadronStringQuality'),
            selection = cms.string("algoIs(\'kPFNeutralHadron\')"),
            selectionFailValue = cms.double(1000.0),
            selectionPassFunction = cms.string('-pt')
        )
    ),
    verbosity = cms.int32(0)
)


process.beamhaloTrackCandidates = cms.EDProducer("CkfTrackCandidateMaker",
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('BeamHaloNavigationSchool'),
    RedundantSeedCleaner = cms.string('CachingSeedCleanerBySharedInput'),
    TrajectoryBuilderPSet = cms.PSet(
        refToPSet_ = cms.string('CkfTrajectoryBuilderBeamHalo')
    ),
    TrajectoryCleaner = cms.string('TrajectoryCleanerBySharedHits'),
    TransientInitialStateEstimatorParameters = cms.PSet(
        numberMeasurementsForFit = cms.int32(4),
        propagatorAlongTISE = cms.string('BeamHaloPropagatorAlong'),
        propagatorOppositeTISE = cms.string('BeamHaloPropagatorOpposite')
    ),
    cleanTrajectoryAfterInOut = cms.bool(True),
    clustersToSkip = cms.InputTag(""),
    doSeedingRegionRebuilding = cms.bool(True),
    maxNSeeds = cms.uint32(500000),
    maxSeedsBeforeCleaning = cms.uint32(5000),
    mightGet = cms.optional.untracked.vstring,
    numHitsForSeedCleaner = cms.int32(4),
    onlyPixelHitsForSeedCleaner = cms.bool(False),
    phase2clustersToSkip = cms.InputTag(""),
    reverseTrajectories = cms.bool(False),
    src = cms.InputTag("beamhaloTrackerSeeds"),
    useHitsSplitting = cms.bool(True)
)


process.beamhaloTrackerSeedingLayers = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(

    ),
    FPix = cms.PSet(
        HitProducer = cms.string('siPixelRecHits'),
        TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4PixelPairs')
    ),
    MTEC = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        maxRing = cms.int32(7),
        minRing = cms.int32(5),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
        useRingSlector = cms.bool(False)
    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        maxRing = cms.int32(2),
        minRing = cms.int32(1),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
        useRingSlector = cms.bool(True)
    ),
    TIB = cms.PSet(

    ),
    TID = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
        useRingSlector = cms.bool(False)
    ),
    TOB = cms.PSet(

    ),
    layerList = cms.vstring(
        'FPix1_pos+FPix2_pos',
        'FPix1_neg+FPix2_neg',
        'TID2_pos+TID3_pos',
        'TID2_neg+TID3_neg',
        'TEC1_neg+TEC2_neg',
        'TEC1_pos+TEC2_pos',
        'TEC2_neg+TEC3_neg',
        'TEC2_pos+TEC3_pos',
        'TEC3_neg+TEC4_neg',
        'TEC3_pos+TEC4_pos',
        'TEC4_neg+TEC5_neg',
        'TEC4_pos+TEC5_pos',
        'TEC5_neg+TEC6_neg',
        'TEC5_pos+TEC6_pos',
        'MTEC7_neg+MTEC8_neg',
        'MTEC7_pos+MTEC8_pos',
        'MTEC8_neg+MTEC9_neg',
        'MTEC8_pos+MTEC9_pos'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.beamhaloTrackerSeeds = cms.EDProducer("CtfSpecialSeedGenerator",
    Charges = cms.vint32(-1, 1),
    CheckHitsAreOnDifferentLayers = cms.bool(False),
    ClusterCollectionLabel = cms.InputTag("siStripClusters"),
    ErrorRescaling = cms.double(50.0),
    MaxNumberOfCosmicClusters = cms.uint32(10000),
    MaxNumberOfPixelClusters = cms.uint32(10000),
    OrderedHitsFactoryPSets = cms.VPSet(
        cms.PSet(
            ComponentName = cms.string('BeamHaloPairGenerator'),
            LayerSrc = cms.InputTag("beamhaloTrackerSeedingLayers"),
            NavigationDirection = cms.string('outsideIn'),
            PropagationDirection = cms.string('alongMomentum'),
            maxTheta = cms.double(0.1)
        ),
        cms.PSet(
            ComponentName = cms.string('BeamHaloPairGenerator'),
            LayerSrc = cms.InputTag("beamhaloTrackerSeedingLayers"),
            NavigationDirection = cms.string('outsideIn'),
            PropagationDirection = cms.string('oppositeToMomentum'),
            maxTheta = cms.double(0.1)
        )
    ),
    PixelClusterCollectionLabel = cms.InputTag("siPixelClusters"),
    RegionFactoryPSet = cms.PSet(
        ComponentName = cms.string('GlobalRegionProducer'),
        RegionPSet = cms.PSet(
            originHalfLength = cms.double(21.2),
            originRadius = cms.double(0.2),
            originXPos = cms.double(0.0),
            originYPos = cms.double(0.0),
            originZPos = cms.double(0.0),
            precise = cms.bool(True),
            ptMin = cms.double(0.9),
            useMultipleScattering = cms.bool(False)
        )
    ),
    SeedMomentum = cms.double(15.0),
    SeedsFromNegativeY = cms.bool(False),
    SeedsFromPositiveY = cms.bool(False),
    SetMomentum = cms.bool(True),
    TTRHBuilder = cms.string('WithTrackAngle'),
    UseScintillatorsConstraint = cms.bool(False),
    doClusterCheck = cms.bool(True),
    maxSeeds = cms.int32(10000),
    requireBOFF = cms.bool(False)
)


process.beamhaloTracks = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('beamhalo'),
    Fitter = cms.string('KFFittingSmootherBH'),
    GeometricInnerState = cms.bool(True),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('BeamHaloNavigationSchool'),
    Propagator = cms.string('BeamHaloPropagatorAlong'),
    SimpleMagneticField = cms.string(''),
    TTRHBuilder = cms.string('WithTrackAngle'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('beamhaloTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("beamhaloTrackCandidates"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(False)
)


process.caloTowerForTrk = cms.EDProducer("CaloTowersCreator",
    AllowMissingInputs = cms.bool(False),
    EBGrid = cms.vdouble(-1.0, 1.0, 10.0, 100.0, 1000.0),
    EBSumThreshold = cms.double(0.2),
    EBThreshold = cms.double(0.07),
    EBWeight = cms.double(1.0),
    EBWeights = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0),
    EEGrid = cms.vdouble(-1.0, 1.0, 10.0, 100.0, 1000.0),
    EESumThreshold = cms.double(0.45),
    EEThreshold = cms.double(0.3),
    EEWeight = cms.double(1.0),
    EEWeights = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0),
    EcalRecHitSeveritiesToBeExcluded = cms.vstring(
        'kTime',
        'kWeird',
        'kBad'
    ),
    EcalSeveritiesToBeUsedInBadTowers = cms.vstring(),
    EcutTower = cms.double(-1000.0),
    HBGrid = cms.vdouble(-1.0, 1.0, 10.0, 100.0, 1000.0),
    HBThreshold = cms.double(0.3),
    HBThreshold1 = cms.double(0.1),
    HBThreshold2 = cms.double(0.2),
    HBWeight = cms.double(1.0),
    HBWeights = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0),
    HEDGrid = cms.vdouble(-1.0, 1.0, 10.0, 100.0, 1000.0),
    HEDThreshold = cms.double(0.2),
    HEDThreshold1 = cms.double(0.1),
    HEDWeight = cms.double(1.0),
    HEDWeights = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0),
    HESGrid = cms.vdouble(-1.0, 1.0, 10.0, 100.0, 1000.0),
    HESThreshold = cms.double(0.2),
    HESThreshold1 = cms.double(0.1),
    HESWeight = cms.double(1.0),
    HESWeights = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0),
    HF1Grid = cms.vdouble(-1.0, 1.0, 10.0, 100.0, 1000.0),
    HF1Threshold = cms.double(0.5),
    HF1Weight = cms.double(1.0),
    HF1Weights = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0),
    HF2Grid = cms.vdouble(-1.0, 1.0, 10.0, 100.0, 1000.0),
    HF2Threshold = cms.double(0.85),
    HF2Weight = cms.double(1.0),
    HF2Weights = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0),
    HOGrid = cms.vdouble(-1.0, 1.0, 10.0, 100.0, 1000.0),
    HOThreshold0 = cms.double(1.1),
    HOThresholdMinus1 = cms.double(3.5),
    HOThresholdMinus2 = cms.double(3.5),
    HOThresholdPlus1 = cms.double(3.5),
    HOThresholdPlus2 = cms.double(3.5),
    HOWeight = cms.double(1.0),
    HOWeights = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0),
    HcalAcceptSeverityLevel = cms.uint32(9),
    HcalAcceptSeverityLevelForRejectedHit = cms.uint32(9999),
    HcalPhase = cms.int32(1),
    HcalThreshold = cms.double(-1000.0),
    MomConstrMethod = cms.int32(1),
    MomEBDepth = cms.double(0.3),
    MomEEDepth = cms.double(0.0),
    MomHBDepth = cms.double(0.2),
    MomHEDepth = cms.double(0.4),
    UseEcalRecoveredHits = cms.bool(False),
    UseEtEBTreshold = cms.bool(False),
    UseEtEETreshold = cms.bool(False),
    UseHO = cms.bool(True),
    UseHcalRecoveredHits = cms.bool(True),
    UseRejectedHitsOnly = cms.bool(False),
    UseRejectedRecoveredEcalHits = cms.bool(False),
    UseRejectedRecoveredHcalHits = cms.bool(True),
    UseSymEBTreshold = cms.bool(True),
    UseSymEETreshold = cms.bool(True),
    ecalInputs = cms.VInputTag(cms.InputTag("ecalRecHit","EcalRecHitsEB"), cms.InputTag("ecalRecHit","EcalRecHitsEE")),
    hbheInput = cms.InputTag("hbhereco"),
    hfInput = cms.InputTag("hfreco"),
    hoInput = cms.InputTag("horeco"),
    missingHcalRescaleFactorForEcal = cms.double(1.0)
)


process.caloTowerForTrkPreSplitting = cms.EDProducer("CaloTowersCreator",
    AllowMissingInputs = cms.bool(False),
    EBGrid = cms.vdouble(-1.0, 1.0, 10.0, 100.0, 1000.0),
    EBSumThreshold = cms.double(0.2),
    EBThreshold = cms.double(0.07),
    EBWeight = cms.double(1.0),
    EBWeights = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0),
    EEGrid = cms.vdouble(-1.0, 1.0, 10.0, 100.0, 1000.0),
    EESumThreshold = cms.double(0.45),
    EEThreshold = cms.double(0.3),
    EEWeight = cms.double(1.0),
    EEWeights = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0),
    EcalRecHitSeveritiesToBeExcluded = cms.vstring(
        'kTime',
        'kWeird',
        'kBad'
    ),
    EcalSeveritiesToBeUsedInBadTowers = cms.vstring(),
    EcutTower = cms.double(-1000.0),
    HBGrid = cms.vdouble(-1.0, 1.0, 10.0, 100.0, 1000.0),
    HBThreshold = cms.double(0.3),
    HBThreshold1 = cms.double(0.1),
    HBThreshold2 = cms.double(0.2),
    HBWeight = cms.double(1.0),
    HBWeights = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0),
    HEDGrid = cms.vdouble(-1.0, 1.0, 10.0, 100.0, 1000.0),
    HEDThreshold = cms.double(0.2),
    HEDThreshold1 = cms.double(0.1),
    HEDWeight = cms.double(1.0),
    HEDWeights = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0),
    HESGrid = cms.vdouble(-1.0, 1.0, 10.0, 100.0, 1000.0),
    HESThreshold = cms.double(0.2),
    HESThreshold1 = cms.double(0.1),
    HESWeight = cms.double(1.0),
    HESWeights = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0),
    HF1Grid = cms.vdouble(-1.0, 1.0, 10.0, 100.0, 1000.0),
    HF1Threshold = cms.double(0.5),
    HF1Weight = cms.double(1.0),
    HF1Weights = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0),
    HF2Grid = cms.vdouble(-1.0, 1.0, 10.0, 100.0, 1000.0),
    HF2Threshold = cms.double(0.85),
    HF2Weight = cms.double(1.0),
    HF2Weights = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0),
    HOGrid = cms.vdouble(-1.0, 1.0, 10.0, 100.0, 1000.0),
    HOThreshold0 = cms.double(1.1),
    HOThresholdMinus1 = cms.double(3.5),
    HOThresholdMinus2 = cms.double(3.5),
    HOThresholdPlus1 = cms.double(3.5),
    HOThresholdPlus2 = cms.double(3.5),
    HOWeight = cms.double(1.0),
    HOWeights = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0),
    HcalAcceptSeverityLevel = cms.uint32(9),
    HcalAcceptSeverityLevelForRejectedHit = cms.uint32(9999),
    HcalPhase = cms.int32(1),
    HcalThreshold = cms.double(-1000.0),
    MomConstrMethod = cms.int32(1),
    MomEBDepth = cms.double(0.3),
    MomEEDepth = cms.double(0.0),
    MomHBDepth = cms.double(0.2),
    MomHEDepth = cms.double(0.4),
    UseEcalRecoveredHits = cms.bool(False),
    UseEtEBTreshold = cms.bool(False),
    UseEtEETreshold = cms.bool(False),
    UseHO = cms.bool(True),
    UseHcalRecoveredHits = cms.bool(True),
    UseRejectedHitsOnly = cms.bool(False),
    UseRejectedRecoveredEcalHits = cms.bool(False),
    UseRejectedRecoveredHcalHits = cms.bool(True),
    UseSymEBTreshold = cms.bool(True),
    UseSymEETreshold = cms.bool(True),
    ecalInputs = cms.VInputTag(cms.InputTag("ecalRecHit","EcalRecHitsEB"), cms.InputTag("ecalRecHit","EcalRecHitsEE")),
    hbheInput = cms.InputTag("hbhereco"),
    hfInput = cms.InputTag("hfreco"),
    hoInput = cms.InputTag("horeco"),
    missingHcalRescaleFactorForEcal = cms.double(1.0)
)


process.calotowermaker = cms.EDProducer("CaloTowersCreator",
    AllowMissingInputs = cms.bool(False),
    EBGrid = cms.vdouble(-1.0, 1.0, 10.0, 100.0, 1000.0),
    EBSumThreshold = cms.double(0.2),
    EBThreshold = cms.double(0.07),
    EBWeight = cms.double(1.0),
    EBWeights = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0),
    EEGrid = cms.vdouble(-1.0, 1.0, 10.0, 100.0, 1000.0),
    EESumThreshold = cms.double(0.45),
    EEThreshold = cms.double(0.3),
    EEWeight = cms.double(1.0),
    EEWeights = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0),
    EcalRecHitSeveritiesToBeExcluded = cms.vstring(
        'kTime',
        'kWeird',
        'kBad'
    ),
    EcalSeveritiesToBeUsedInBadTowers = cms.vstring(),
    EcutTower = cms.double(-1000.0),
    HBGrid = cms.vdouble(-1.0, 1.0, 10.0, 100.0, 1000.0),
    HBThreshold = cms.double(0.3),
    HBThreshold1 = cms.double(0.1),
    HBThreshold2 = cms.double(0.2),
    HBWeight = cms.double(1.0),
    HBWeights = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0),
    HEDGrid = cms.vdouble(-1.0, 1.0, 10.0, 100.0, 1000.0),
    HEDThreshold = cms.double(0.2),
    HEDThreshold1 = cms.double(0.1),
    HEDWeight = cms.double(1.0),
    HEDWeights = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0),
    HESGrid = cms.vdouble(-1.0, 1.0, 10.0, 100.0, 1000.0),
    HESThreshold = cms.double(0.2),
    HESThreshold1 = cms.double(0.1),
    HESWeight = cms.double(1.0),
    HESWeights = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0),
    HF1Grid = cms.vdouble(-1.0, 1.0, 10.0, 100.0, 1000.0),
    HF1Threshold = cms.double(0.5),
    HF1Weight = cms.double(1.0),
    HF1Weights = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0),
    HF2Grid = cms.vdouble(-1.0, 1.0, 10.0, 100.0, 1000.0),
    HF2Threshold = cms.double(0.85),
    HF2Weight = cms.double(1.0),
    HF2Weights = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0),
    HOGrid = cms.vdouble(-1.0, 1.0, 10.0, 100.0, 1000.0),
    HOThreshold0 = cms.double(1.1),
    HOThresholdMinus1 = cms.double(3.5),
    HOThresholdMinus2 = cms.double(3.5),
    HOThresholdPlus1 = cms.double(3.5),
    HOThresholdPlus2 = cms.double(3.5),
    HOWeight = cms.double(1.0),
    HOWeights = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0),
    HcalAcceptSeverityLevel = cms.uint32(9),
    HcalAcceptSeverityLevelForRejectedHit = cms.uint32(9999),
    HcalPhase = cms.int32(1),
    HcalThreshold = cms.double(-1000.0),
    MomConstrMethod = cms.int32(1),
    MomEBDepth = cms.double(0.3),
    MomEEDepth = cms.double(0.0),
    MomHBDepth = cms.double(0.2),
    MomHEDepth = cms.double(0.4),
    UseEcalRecoveredHits = cms.bool(False),
    UseEtEBTreshold = cms.bool(False),
    UseEtEETreshold = cms.bool(False),
    UseHO = cms.bool(True),
    UseHcalRecoveredHits = cms.bool(True),
    UseRejectedHitsOnly = cms.bool(False),
    UseRejectedRecoveredEcalHits = cms.bool(False),
    UseRejectedRecoveredHcalHits = cms.bool(True),
    UseSymEBTreshold = cms.bool(True),
    UseSymEETreshold = cms.bool(True),
    ecalInputs = cms.VInputTag(cms.InputTag("ecalRecHit","EcalRecHitsEB"), cms.InputTag("ecalRecHit","EcalRecHitsEE")),
    hbheInput = cms.InputTag("hbhereco"),
    hfInput = cms.InputTag("hfreco"),
    hoInput = cms.InputTag("horeco"),
    missingHcalRescaleFactorForEcal = cms.double(0)
)


process.chargeCut2069Clusters = cms.EDProducer("ClusterChargeMasker",
    clusterChargeCut = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutTight')
    ),
    oldClusterRemovalInfo = cms.InputTag("pixelPairStepClusters"),
    pixelClusters = cms.InputTag("siPixelClusters"),
    stripClusters = cms.InputTag("siStripClusters")
)


process.ckfTrackCandidates = cms.EDProducer("CkfTrackCandidateMaker",
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    RedundantSeedCleaner = cms.string('CachingSeedCleanerBySharedInput'),
    TrajectoryBuilderPSet = cms.PSet(
        refToPSet_ = cms.string('GroupedCkfTrajectoryBuilder')
    ),
    TrajectoryCleaner = cms.string('TrajectoryCleanerBySharedHits'),
    TransientInitialStateEstimatorParameters = cms.PSet(
        numberMeasurementsForFit = cms.int32(4),
        propagatorAlongTISE = cms.string('PropagatorWithMaterial'),
        propagatorOppositeTISE = cms.string('PropagatorWithMaterialOpposite')
    ),
    cleanTrajectoryAfterInOut = cms.bool(True),
    clustersToSkip = cms.InputTag(""),
    doSeedingRegionRebuilding = cms.bool(True),
    maxNSeeds = cms.uint32(500000),
    maxSeedsBeforeCleaning = cms.uint32(5000),
    mightGet = cms.optional.untracked.vstring,
    numHitsForSeedCleaner = cms.int32(4),
    onlyPixelHitsForSeedCleaner = cms.bool(False),
    phase2clustersToSkip = cms.InputTag(""),
    reverseTrajectories = cms.bool(False),
    src = cms.InputTag("globalMixedSeeds"),
    useHitsSplitting = cms.bool(True)
)


process.ckfTrackCandidatesCombinedSeeds = cms.EDProducer("CkfTrackCandidateMaker",
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    RedundantSeedCleaner = cms.string('CachingSeedCleanerBySharedInput'),
    TrajectoryBuilderPSet = cms.PSet(
        refToPSet_ = cms.string('GroupedCkfTrajectoryBuilder')
    ),
    TrajectoryCleaner = cms.string('TrajectoryCleanerBySharedHits'),
    TransientInitialStateEstimatorParameters = cms.PSet(
        numberMeasurementsForFit = cms.int32(4),
        propagatorAlongTISE = cms.string('PropagatorWithMaterial'),
        propagatorOppositeTISE = cms.string('PropagatorWithMaterialOpposite')
    ),
    cleanTrajectoryAfterInOut = cms.bool(True),
    clustersToSkip = cms.InputTag(""),
    doSeedingRegionRebuilding = cms.bool(True),
    maxNSeeds = cms.uint32(500000),
    maxSeedsBeforeCleaning = cms.uint32(5000),
    mightGet = cms.optional.untracked.vstring,
    numHitsForSeedCleaner = cms.int32(4),
    onlyPixelHitsForSeedCleaner = cms.bool(False),
    phase2clustersToSkip = cms.InputTag(""),
    reverseTrajectories = cms.bool(False),
    src = cms.InputTag("globalCombinedSeeds"),
    useHitsSplitting = cms.bool(True)
)


process.ckfTrackCandidatesNoOverlaps = cms.EDProducer("CkfTrackCandidateMaker",
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    RedundantSeedCleaner = cms.string('CachingSeedCleanerBySharedInput'),
    TrajectoryBuilderPSet = cms.PSet(
        refToPSet_ = cms.string('CkfTrajectoryBuilder')
    ),
    TrajectoryCleaner = cms.string('TrajectoryCleanerBySharedHits'),
    TransientInitialStateEstimatorParameters = cms.PSet(
        numberMeasurementsForFit = cms.int32(4),
        propagatorAlongTISE = cms.string('PropagatorWithMaterial'),
        propagatorOppositeTISE = cms.string('PropagatorWithMaterialOpposite')
    ),
    cleanTrajectoryAfterInOut = cms.bool(True),
    clustersToSkip = cms.InputTag(""),
    doSeedingRegionRebuilding = cms.bool(True),
    maxNSeeds = cms.uint32(500000),
    maxSeedsBeforeCleaning = cms.uint32(5000),
    mightGet = cms.optional.untracked.vstring,
    numHitsForSeedCleaner = cms.int32(4),
    onlyPixelHitsForSeedCleaner = cms.bool(False),
    phase2clustersToSkip = cms.InputTag(""),
    reverseTrajectories = cms.bool(False),
    src = cms.InputTag("globalMixedSeeds"),
    useHitsSplitting = cms.bool(True)
)


process.ckfTrackCandidatesP5 = cms.EDProducer("CkfTrackCandidateMaker",
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('CosmicNavigationSchool'),
    RedundantSeedCleaner = cms.string('CachingSeedCleanerBySharedInput'),
    TrajectoryBuilderPSet = cms.PSet(
        refToPSet_ = cms.string('GroupedCkfTrajectoryBuilderP5')
    ),
    TrajectoryCleaner = cms.string('TrajectoryCleanerBySharedHits'),
    TransientInitialStateEstimatorParameters = cms.PSet(
        numberMeasurementsForFit = cms.int32(4),
        propagatorAlongTISE = cms.string('PropagatorWithMaterial'),
        propagatorOppositeTISE = cms.string('PropagatorWithMaterialOpposite')
    ),
    cleanTrajectoryAfterInOut = cms.bool(True),
    clustersToSkip = cms.InputTag(""),
    doSeedingRegionRebuilding = cms.bool(True),
    maxNSeeds = cms.uint32(500000),
    maxSeedsBeforeCleaning = cms.uint32(5000),
    mightGet = cms.optional.untracked.vstring,
    numHitsForSeedCleaner = cms.int32(4),
    onlyPixelHitsForSeedCleaner = cms.bool(False),
    phase2clustersToSkip = cms.InputTag(""),
    reverseTrajectories = cms.bool(False),
    src = cms.InputTag("combinedP5SeedsForCTF"),
    useHitsSplitting = cms.bool(True)
)


process.ckfTrackCandidatesP5Bottom = cms.EDProducer("CkfTrackCandidateMaker",
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('CosmicNavigationSchool'),
    RedundantSeedCleaner = cms.string('CachingSeedCleanerBySharedInput'),
    TrajectoryBuilderPSet = cms.PSet(
        refToPSet_ = cms.string('GroupedCkfTrajectoryBuilderP5Bottom')
    ),
    TrajectoryCleaner = cms.string('TrajectoryCleanerBySharedHits'),
    TransientInitialStateEstimatorParameters = cms.PSet(
        numberMeasurementsForFit = cms.int32(4),
        propagatorAlongTISE = cms.string('PropagatorWithMaterial'),
        propagatorOppositeTISE = cms.string('PropagatorWithMaterialOpposite')
    ),
    cleanTrajectoryAfterInOut = cms.bool(True),
    clustersToSkip = cms.InputTag(""),
    doSeedingRegionRebuilding = cms.bool(True),
    maxNSeeds = cms.uint32(500000),
    maxSeedsBeforeCleaning = cms.uint32(5000),
    mightGet = cms.optional.untracked.vstring,
    numHitsForSeedCleaner = cms.int32(4),
    onlyPixelHitsForSeedCleaner = cms.bool(False),
    phase2clustersToSkip = cms.InputTag(""),
    reverseTrajectories = cms.bool(False),
    src = cms.InputTag("combinedP5SeedsForCTFBottom"),
    useHitsSplitting = cms.bool(True)
)


process.ckfTrackCandidatesP5LHCNavigation = cms.EDProducer("CkfTrackCandidateMaker",
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    RedundantSeedCleaner = cms.string('CachingSeedCleanerBySharedInput'),
    TrajectoryBuilderPSet = cms.PSet(
        refToPSet_ = cms.string('GroupedCkfTrajectoryBuilderP5')
    ),
    TrajectoryCleaner = cms.string('TrajectoryCleanerBySharedHits'),
    TransientInitialStateEstimatorParameters = cms.PSet(
        numberMeasurementsForFit = cms.int32(4),
        propagatorAlongTISE = cms.string('PropagatorWithMaterial'),
        propagatorOppositeTISE = cms.string('PropagatorWithMaterialOpposite')
    ),
    cleanTrajectoryAfterInOut = cms.bool(True),
    clustersToSkip = cms.InputTag(""),
    doSeedingRegionRebuilding = cms.bool(True),
    maxNSeeds = cms.uint32(500000),
    maxSeedsBeforeCleaning = cms.uint32(5000),
    mightGet = cms.optional.untracked.vstring,
    numHitsForSeedCleaner = cms.int32(4),
    onlyPixelHitsForSeedCleaner = cms.bool(False),
    phase2clustersToSkip = cms.InputTag(""),
    reverseTrajectories = cms.bool(False),
    src = cms.InputTag("combinedP5SeedsForCTF"),
    useHitsSplitting = cms.bool(True)
)


process.ckfTrackCandidatesP5Top = cms.EDProducer("CkfTrackCandidateMaker",
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('CosmicNavigationSchool'),
    RedundantSeedCleaner = cms.string('CachingSeedCleanerBySharedInput'),
    TrajectoryBuilderPSet = cms.PSet(
        refToPSet_ = cms.string('GroupedCkfTrajectoryBuilderP5Top')
    ),
    TrajectoryCleaner = cms.string('TrajectoryCleanerBySharedHits'),
    TransientInitialStateEstimatorParameters = cms.PSet(
        numberMeasurementsForFit = cms.int32(4),
        propagatorAlongTISE = cms.string('PropagatorWithMaterial'),
        propagatorOppositeTISE = cms.string('PropagatorWithMaterialOpposite')
    ),
    cleanTrajectoryAfterInOut = cms.bool(True),
    clustersToSkip = cms.InputTag(""),
    doSeedingRegionRebuilding = cms.bool(True),
    maxNSeeds = cms.uint32(500000),
    maxSeedsBeforeCleaning = cms.uint32(5000),
    mightGet = cms.optional.untracked.vstring,
    numHitsForSeedCleaner = cms.int32(4),
    onlyPixelHitsForSeedCleaner = cms.bool(False),
    phase2clustersToSkip = cms.InputTag(""),
    reverseTrajectories = cms.bool(False),
    src = cms.InputTag("combinedP5SeedsForCTFTop"),
    useHitsSplitting = cms.bool(True)
)


process.ckfTrackCandidatesPixelLess = cms.EDProducer("CkfTrackCandidateMaker",
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    RedundantSeedCleaner = cms.string('CachingSeedCleanerBySharedInput'),
    TrajectoryBuilderPSet = cms.PSet(
        refToPSet_ = cms.string('GroupedCkfTrajectoryBuilder')
    ),
    TrajectoryCleaner = cms.string('TrajectoryCleanerBySharedHits'),
    TransientInitialStateEstimatorParameters = cms.PSet(
        numberMeasurementsForFit = cms.int32(4),
        propagatorAlongTISE = cms.string('PropagatorWithMaterial'),
        propagatorOppositeTISE = cms.string('PropagatorWithMaterialOpposite')
    ),
    cleanTrajectoryAfterInOut = cms.bool(True),
    clustersToSkip = cms.InputTag(""),
    doSeedingRegionRebuilding = cms.bool(True),
    maxNSeeds = cms.uint32(500000),
    maxSeedsBeforeCleaning = cms.uint32(5000),
    mightGet = cms.optional.untracked.vstring,
    numHitsForSeedCleaner = cms.int32(4),
    onlyPixelHitsForSeedCleaner = cms.bool(False),
    phase2clustersToSkip = cms.InputTag(""),
    reverseTrajectories = cms.bool(False),
    src = cms.InputTag("globalPixelLessSeeds"),
    useHitsSplitting = cms.bool(True)
)


process.clusterSummaryProducer = cms.EDProducer("ClusterSummaryProducer",
    doPixels = cms.bool(True),
    doStrips = cms.bool(True),
    pixelClusters = cms.InputTag("siPixelClustersPreSplitting"),
    stripClusters = cms.InputTag("siStripClusters"),
    verbose = cms.bool(False),
    wantedSubDets = cms.vstring(
        'TOB',
        'TIB',
        'TID',
        'TEC',
        'STRIP',
        'BPIX',
        'FPIX',
        'PIXEL'
    ),
    wantedUserSubDets = cms.VPSet()
)


process.clusterSummaryProducerNoSplitting = cms.EDProducer("ClusterSummaryProducer",
    doPixels = cms.bool(True),
    doStrips = cms.bool(True),
    pixelClusters = cms.InputTag("siPixelClusters"),
    stripClusters = cms.InputTag("siStripClusters"),
    verbose = cms.bool(False),
    wantedSubDets = cms.vstring(
        'TOB',
        'TIB',
        'TID',
        'TEC',
        'STRIP',
        'BPIX',
        'FPIX',
        'PIXEL'
    ),
    wantedUserSubDets = cms.VPSet()
)


process.combinatorialcosmicseedfinder = cms.EDProducer("CtfSpecialSeedGenerator",
    Charges = cms.vint32(-1),
    CheckHitsAreOnDifferentLayers = cms.bool(False),
    ClusterCollectionLabel = cms.InputTag("siStripClusters"),
    DontCountDetsAboveNClusters = cms.uint32(20),
    ErrorRescaling = cms.double(50.0),
    LowerScintillatorParameters = cms.PSet(
        GlobalX = cms.double(0.0),
        GlobalY = cms.double(-100.0),
        GlobalZ = cms.double(50.0),
        LenghtInZ = cms.double(100.0),
        WidthInX = cms.double(100.0)
    ),
    MaxNumberOfCosmicClusters = cms.uint32(300),
    MaxNumberOfPixelClusters = cms.uint32(300),
    OrderedHitsFactoryPSets = cms.VPSet(
        cms.PSet(
            ComponentName = cms.string('GenericTripletGenerator'),
            LayerSrc = cms.InputTag("combinatorialcosmicseedingtripletsTOB"),
            NavigationDirection = cms.string('outsideIn'),
            PropagationDirection = cms.string('alongMomentum')
        ),
        cms.PSet(
            ComponentName = cms.string('GenericPairGenerator'),
            LayerSrc = cms.InputTag("combinatorialcosmicseedingpairsTECpos"),
            NavigationDirection = cms.string('outsideIn'),
            PropagationDirection = cms.string('alongMomentum')
        ),
        cms.PSet(
            ComponentName = cms.string('GenericTripletGenerator'),
            LayerSrc = cms.InputTag("combinatorialcosmicseedingtripletsTIB"),
            NavigationDirection = cms.string('insideOut'),
            PropagationDirection = cms.string('oppositeToMomentum')
        )
    ),
    PixelClusterCollectionLabel = cms.InputTag("siPixelClusters"),
    RegionFactoryPSet = cms.PSet(
        ComponentName = cms.string('GlobalRegionProducer'),
        RegionPSet = cms.PSet(
            originHalfLength = cms.double(21.2),
            originRadius = cms.double(0.2),
            originXPos = cms.double(0.0),
            originYPos = cms.double(0.0),
            originZPos = cms.double(0.0),
            precise = cms.bool(True),
            ptMin = cms.double(0.9),
            useMultipleScattering = cms.bool(False)
        )
    ),
    SeedMomentum = cms.double(5.0),
    SeedsFromNegativeY = cms.bool(False),
    SeedsFromPositiveY = cms.bool(True),
    SetMomentum = cms.bool(True),
    TTRHBuilder = cms.string('WithTrackAngle'),
    UpperScintillatorParameters = cms.PSet(
        GlobalX = cms.double(0.0),
        GlobalY = cms.double(300.0),
        GlobalZ = cms.double(50.0),
        LenghtInZ = cms.double(100.0),
        WidthInX = cms.double(100.0)
    ),
    UseScintillatorsConstraint = cms.bool(True),
    doClusterCheck = cms.bool(True),
    maxSeeds = cms.int32(10000),
    requireBOFF = cms.bool(False)
)


process.combinatorialcosmicseedfinderP5 = cms.EDProducer("CtfSpecialSeedGenerator",
    Charges = cms.vint32(-1),
    CheckHitsAreOnDifferentLayers = cms.bool(False),
    ClusterCollectionLabel = cms.InputTag("siStripClusters"),
    DontCountDetsAboveNClusters = cms.uint32(20),
    ErrorRescaling = cms.double(50.0),
    LowerScintillatorParameters = cms.PSet(
        GlobalX = cms.double(0.0),
        GlobalY = cms.double(-100.0),
        GlobalZ = cms.double(50.0),
        LenghtInZ = cms.double(100.0),
        WidthInX = cms.double(100.0)
    ),
    MaxNumberOfCosmicClusters = cms.uint32(300),
    MaxNumberOfPixelClusters = cms.uint32(300),
    OrderedHitsFactoryPSets = cms.VPSet(
        cms.PSet(
            ComponentName = cms.string('GenericTripletGenerator'),
            LayerSrc = cms.InputTag("combinatorialcosmicseedingtripletsP5"),
            NavigationDirection = cms.string('outsideIn'),
            PropagationDirection = cms.string('alongMomentum')
        ),
        cms.PSet(
            ComponentName = cms.string('GenericPairGenerator'),
            LayerSrc = cms.InputTag("combinatorialcosmicseedingpairsTOBP5"),
            NavigationDirection = cms.string('outsideIn'),
            PropagationDirection = cms.string('alongMomentum')
        ),
        cms.PSet(
            ComponentName = cms.string('GenericPairGenerator'),
            LayerSrc = cms.InputTag("combinatorialcosmicseedingpairsTECposP5"),
            NavigationDirection = cms.string('outsideIn'),
            PropagationDirection = cms.string('alongMomentum')
        ),
        cms.PSet(
            ComponentName = cms.string('GenericPairGenerator'),
            LayerSrc = cms.InputTag("combinatorialcosmicseedingpairsTECposP5"),
            NavigationDirection = cms.string('insideOut'),
            PropagationDirection = cms.string('alongMomentum')
        ),
        cms.PSet(
            ComponentName = cms.string('GenericPairGenerator'),
            LayerSrc = cms.InputTag("combinatorialcosmicseedingpairsTECnegP5"),
            NavigationDirection = cms.string('outsideIn'),
            PropagationDirection = cms.string('alongMomentum')
        ),
        cms.PSet(
            ComponentName = cms.string('GenericPairGenerator'),
            LayerSrc = cms.InputTag("combinatorialcosmicseedingpairsTECnegP5"),
            NavigationDirection = cms.string('insideOut'),
            PropagationDirection = cms.string('alongMomentum')
        )
    ),
    PixelClusterCollectionLabel = cms.InputTag("siPixelClusters"),
    RegionFactoryPSet = cms.PSet(
        ComponentName = cms.string('GlobalRegionProducer'),
        RegionPSet = cms.PSet(
            originHalfLength = cms.double(21.2),
            originRadius = cms.double(0.2),
            originXPos = cms.double(0.0),
            originYPos = cms.double(0.0),
            originZPos = cms.double(0.0),
            precise = cms.bool(True),
            ptMin = cms.double(0.9),
            useMultipleScattering = cms.bool(False)
        )
    ),
    SeedMomentum = cms.double(5.0),
    SeedsFromNegativeY = cms.bool(False),
    SeedsFromPositiveY = cms.bool(True),
    SetMomentum = cms.bool(True),
    TTRHBuilder = cms.string('WithTrackAngle'),
    UpperScintillatorParameters = cms.PSet(
        GlobalX = cms.double(0.0),
        GlobalY = cms.double(300.0),
        GlobalZ = cms.double(50.0),
        LenghtInZ = cms.double(100.0),
        WidthInX = cms.double(100.0)
    ),
    UseScintillatorsConstraint = cms.bool(False),
    doClusterCheck = cms.bool(True),
    maxSeeds = cms.int32(10000),
    requireBOFF = cms.bool(True)
)


process.combinatorialcosmicseedfinderP5Bottom = cms.EDProducer("CtfSpecialSeedGenerator",
    Charges = cms.vint32(-1),
    CheckHitsAreOnDifferentLayers = cms.bool(False),
    ClusterCollectionLabel = cms.InputTag("siStripClustersBottom"),
    DontCountDetsAboveNClusters = cms.uint32(20),
    ErrorRescaling = cms.double(50.0),
    LowerScintillatorParameters = cms.PSet(
        GlobalX = cms.double(0.0),
        GlobalY = cms.double(-100.0),
        GlobalZ = cms.double(50.0),
        LenghtInZ = cms.double(100.0),
        WidthInX = cms.double(100.0)
    ),
    MaxNumberOfCosmicClusters = cms.uint32(150),
    MaxNumberOfPixelClusters = cms.uint32(300),
    OrderedHitsFactoryPSets = cms.VPSet(
        cms.PSet(
            ComponentName = cms.string('GenericTripletGenerator'),
            LayerSrc = cms.InputTag("combinatorialcosmicseedingtripletsP5Bottom"),
            NavigationDirection = cms.string('outsideIn'),
            PropagationDirection = cms.string('oppositeToMomentum')
        ),
        cms.PSet(
            ComponentName = cms.string('GenericPairGenerator'),
            LayerSrc = cms.InputTag("combinatorialcosmicseedingpairsTOBP5Bottom"),
            NavigationDirection = cms.string('outsideIn'),
            PropagationDirection = cms.string('oppositeToMomentum')
        ),
        cms.PSet(
            ComponentName = cms.string('GenericPairGenerator'),
            LayerSrc = cms.InputTag("combinatorialcosmicseedingpairsTECposP5Bottom"),
            NavigationDirection = cms.string('outsideIn'),
            PropagationDirection = cms.string('oppositeToMomentum')
        ),
        cms.PSet(
            ComponentName = cms.string('GenericPairGenerator'),
            LayerSrc = cms.InputTag("combinatorialcosmicseedingpairsTECposP5Bottom"),
            NavigationDirection = cms.string('insideOut'),
            PropagationDirection = cms.string('oppositeToMomentum')
        ),
        cms.PSet(
            ComponentName = cms.string('GenericPairGenerator'),
            LayerSrc = cms.InputTag("combinatorialcosmicseedingpairsTECnegP5Bottom"),
            NavigationDirection = cms.string('outsideIn'),
            PropagationDirection = cms.string('oppositeToMomentum')
        ),
        cms.PSet(
            ComponentName = cms.string('GenericPairGenerator'),
            LayerSrc = cms.InputTag("combinatorialcosmicseedingpairsTECnegP5Bottom"),
            NavigationDirection = cms.string('insideOut'),
            PropagationDirection = cms.string('oppositeToMomentum')
        )
    ),
    PixelClusterCollectionLabel = cms.InputTag("siPixelClusters"),
    RegionFactoryPSet = cms.PSet(
        ComponentName = cms.string('GlobalRegionProducer'),
        RegionPSet = cms.PSet(
            originHalfLength = cms.double(21.2),
            originRadius = cms.double(0.2),
            originXPos = cms.double(0.0),
            originYPos = cms.double(0.0),
            originZPos = cms.double(0.0),
            precise = cms.bool(True),
            ptMin = cms.double(0.9),
            useMultipleScattering = cms.bool(False)
        )
    ),
    SeedMomentum = cms.double(5.0),
    SeedsFromNegativeY = cms.bool(True),
    SeedsFromPositiveY = cms.bool(False),
    SetMomentum = cms.bool(True),
    TTRHBuilder = cms.string('WithTrackAngle'),
    UpperScintillatorParameters = cms.PSet(
        GlobalX = cms.double(0.0),
        GlobalY = cms.double(300.0),
        GlobalZ = cms.double(50.0),
        LenghtInZ = cms.double(100.0),
        WidthInX = cms.double(100.0)
    ),
    UseScintillatorsConstraint = cms.bool(False),
    doClusterCheck = cms.bool(True),
    maxSeeds = cms.int32(10000),
    requireBOFF = cms.bool(True)
)


process.combinatorialcosmicseedfinderP5Top = cms.EDProducer("CtfSpecialSeedGenerator",
    Charges = cms.vint32(-1),
    CheckHitsAreOnDifferentLayers = cms.bool(False),
    ClusterCollectionLabel = cms.InputTag("siStripClustersTop"),
    DontCountDetsAboveNClusters = cms.uint32(20),
    ErrorRescaling = cms.double(50.0),
    LowerScintillatorParameters = cms.PSet(
        GlobalX = cms.double(0.0),
        GlobalY = cms.double(-100.0),
        GlobalZ = cms.double(50.0),
        LenghtInZ = cms.double(100.0),
        WidthInX = cms.double(100.0)
    ),
    MaxNumberOfCosmicClusters = cms.uint32(150),
    MaxNumberOfPixelClusters = cms.uint32(300),
    OrderedHitsFactoryPSets = cms.VPSet(
        cms.PSet(
            ComponentName = cms.string('GenericTripletGenerator'),
            LayerSrc = cms.InputTag("combinatorialcosmicseedingtripletsP5Top"),
            NavigationDirection = cms.string('outsideIn'),
            PropagationDirection = cms.string('alongMomentum')
        ),
        cms.PSet(
            ComponentName = cms.string('GenericPairGenerator'),
            LayerSrc = cms.InputTag("combinatorialcosmicseedingpairsTOBP5Top"),
            NavigationDirection = cms.string('outsideIn'),
            PropagationDirection = cms.string('alongMomentum')
        ),
        cms.PSet(
            ComponentName = cms.string('GenericPairGenerator'),
            LayerSrc = cms.InputTag("combinatorialcosmicseedingpairsTECposP5Top"),
            NavigationDirection = cms.string('outsideIn'),
            PropagationDirection = cms.string('alongMomentum')
        ),
        cms.PSet(
            ComponentName = cms.string('GenericPairGenerator'),
            LayerSrc = cms.InputTag("combinatorialcosmicseedingpairsTECposP5Top"),
            NavigationDirection = cms.string('insideOut'),
            PropagationDirection = cms.string('alongMomentum')
        ),
        cms.PSet(
            ComponentName = cms.string('GenericPairGenerator'),
            LayerSrc = cms.InputTag("combinatorialcosmicseedingpairsTECnegP5Top"),
            NavigationDirection = cms.string('outsideIn'),
            PropagationDirection = cms.string('alongMomentum')
        ),
        cms.PSet(
            ComponentName = cms.string('GenericPairGenerator'),
            LayerSrc = cms.InputTag("combinatorialcosmicseedingpairsTECnegP5Top"),
            NavigationDirection = cms.string('insideOut'),
            PropagationDirection = cms.string('alongMomentum')
        )
    ),
    PixelClusterCollectionLabel = cms.InputTag("siPixelClusters"),
    RegionFactoryPSet = cms.PSet(
        ComponentName = cms.string('GlobalRegionProducer'),
        RegionPSet = cms.PSet(
            originHalfLength = cms.double(21.2),
            originRadius = cms.double(0.2),
            originXPos = cms.double(0.0),
            originYPos = cms.double(0.0),
            originZPos = cms.double(0.0),
            precise = cms.bool(True),
            ptMin = cms.double(0.9),
            useMultipleScattering = cms.bool(False)
        )
    ),
    SeedMomentum = cms.double(5.0),
    SeedsFromNegativeY = cms.bool(False),
    SeedsFromPositiveY = cms.bool(True),
    SetMomentum = cms.bool(True),
    TTRHBuilder = cms.string('WithTrackAngle'),
    UpperScintillatorParameters = cms.PSet(
        GlobalX = cms.double(0.0),
        GlobalY = cms.double(300.0),
        GlobalZ = cms.double(50.0),
        LenghtInZ = cms.double(100.0),
        WidthInX = cms.double(100.0)
    ),
    UseScintillatorsConstraint = cms.bool(False),
    doClusterCheck = cms.bool(True),
    maxSeeds = cms.int32(10000),
    requireBOFF = cms.bool(True)
)


process.combinatorialcosmicseedingpairsTECnegP5 = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(

    ),
    FPix = cms.PSet(

    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        maxRing = cms.int32(7),
        minRing = cms.int32(5),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
        useRingSlector = cms.bool(True)
    ),
    TIB = cms.PSet(

    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(

    ),
    layerList = cms.vstring(
        'TEC1_neg+TEC2_neg',
        'TEC2_neg+TEC3_neg',
        'TEC3_neg+TEC4_neg',
        'TEC4_neg+TEC5_neg',
        'TEC5_neg+TEC6_neg',
        'TEC6_neg+TEC7_neg',
        'TEC7_neg+TEC8_neg',
        'TEC8_neg+TEC9_neg'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.combinatorialcosmicseedingpairsTECnegP5Bottom = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(

    ),
    FPix = cms.PSet(

    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHitsBottom","matchedRecHit"),
        maxRing = cms.int32(7),
        minRing = cms.int32(5),
        rphiRecHits = cms.InputTag("siStripMatchedRecHitsBottom","rphiRecHit"),
        useRingSlector = cms.bool(True)
    ),
    TIB = cms.PSet(

    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(

    ),
    layerList = cms.vstring(
        'TEC1_neg+TEC2_neg',
        'TEC2_neg+TEC3_neg',
        'TEC3_neg+TEC4_neg',
        'TEC4_neg+TEC5_neg',
        'TEC5_neg+TEC6_neg',
        'TEC6_neg+TEC7_neg',
        'TEC7_neg+TEC8_neg',
        'TEC8_neg+TEC9_neg'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.combinatorialcosmicseedingpairsTECnegP5Top = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(

    ),
    FPix = cms.PSet(

    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHitsTop","matchedRecHit"),
        maxRing = cms.int32(7),
        minRing = cms.int32(5),
        rphiRecHits = cms.InputTag("siStripMatchedRecHitsTop","rphiRecHit"),
        useRingSlector = cms.bool(True)
    ),
    TIB = cms.PSet(

    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(

    ),
    layerList = cms.vstring(
        'TEC1_neg+TEC2_neg',
        'TEC2_neg+TEC3_neg',
        'TEC3_neg+TEC4_neg',
        'TEC4_neg+TEC5_neg',
        'TEC5_neg+TEC6_neg',
        'TEC6_neg+TEC7_neg',
        'TEC7_neg+TEC8_neg',
        'TEC8_neg+TEC9_neg'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.combinatorialcosmicseedingpairsTECposP5 = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(

    ),
    FPix = cms.PSet(

    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        maxRing = cms.int32(7),
        minRing = cms.int32(5),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
        useRingSlector = cms.bool(True)
    ),
    TIB = cms.PSet(

    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(

    ),
    layerList = cms.vstring(
        'TEC1_pos+TEC2_pos',
        'TEC2_pos+TEC3_pos',
        'TEC3_pos+TEC4_pos',
        'TEC4_pos+TEC5_pos',
        'TEC5_pos+TEC6_pos',
        'TEC6_pos+TEC7_pos',
        'TEC7_pos+TEC8_pos',
        'TEC8_pos+TEC9_pos'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.combinatorialcosmicseedingpairsTECposP5Bottom = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(

    ),
    FPix = cms.PSet(

    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHitsBottom","matchedRecHit"),
        maxRing = cms.int32(7),
        minRing = cms.int32(5),
        rphiRecHits = cms.InputTag("siStripMatchedRecHitsBottom","rphiRecHit"),
        useRingSlector = cms.bool(True)
    ),
    TIB = cms.PSet(

    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(

    ),
    layerList = cms.vstring(
        'TEC1_pos+TEC2_pos',
        'TEC2_pos+TEC3_pos',
        'TEC3_pos+TEC4_pos',
        'TEC4_pos+TEC5_pos',
        'TEC5_pos+TEC6_pos',
        'TEC6_pos+TEC7_pos',
        'TEC7_pos+TEC8_pos',
        'TEC8_pos+TEC9_pos'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.combinatorialcosmicseedingpairsTECposP5Top = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(

    ),
    FPix = cms.PSet(

    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHitsTop","matchedRecHit"),
        maxRing = cms.int32(7),
        minRing = cms.int32(5),
        rphiRecHits = cms.InputTag("siStripMatchedRecHitsTop","rphiRecHit"),
        useRingSlector = cms.bool(True)
    ),
    TIB = cms.PSet(

    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(

    ),
    layerList = cms.vstring(
        'TEC1_pos+TEC2_pos',
        'TEC2_pos+TEC3_pos',
        'TEC3_pos+TEC4_pos',
        'TEC4_pos+TEC5_pos',
        'TEC5_pos+TEC6_pos',
        'TEC6_pos+TEC7_pos',
        'TEC7_pos+TEC8_pos',
        'TEC8_pos+TEC9_pos'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.combinatorialcosmicseedingpairsTOBP5 = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(

    ),
    FPix = cms.PSet(

    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit")
    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit")
    ),
    TEC = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        maxRing = cms.int32(7),
        minRing = cms.int32(5),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
        useRingSlector = cms.bool(False),
        useSimpleRphiHitsCleaner = cms.bool(True)
    ),
    TIB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit")
    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit")
    ),
    layerList = cms.vstring(
        'MTOB5+MTOB6',
        'MTOB4+MTOB5'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.combinatorialcosmicseedingpairsTOBP5Bottom = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(

    ),
    FPix = cms.PSet(

    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        rphiRecHits = cms.InputTag("siStripMatchedRecHitsBottom","rphiRecHit")
    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        rphiRecHits = cms.InputTag("siStripMatchedRecHitsBottom","rphiRecHit")
    ),
    TEC = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHitsBottom","matchedRecHit"),
        maxRing = cms.int32(7),
        minRing = cms.int32(5),
        rphiRecHits = cms.InputTag("siStripMatchedRecHitsBottom","rphiRecHit"),
        useRingSlector = cms.bool(False),
        useSimpleRphiHitsCleaner = cms.bool(True)
    ),
    TIB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHitsBottom","matchedRecHit")
    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHitsBottom","matchedRecHit")
    ),
    layerList = cms.vstring(
        'MTOB5+MTOB6',
        'MTOB4+MTOB5'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.combinatorialcosmicseedingpairsTOBP5Top = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(

    ),
    FPix = cms.PSet(

    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        rphiRecHits = cms.InputTag("siStripMatchedRecHitsTop","rphiRecHit")
    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        rphiRecHits = cms.InputTag("siStripMatchedRecHitsTop","rphiRecHit")
    ),
    TEC = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHitsTop","matchedRecHit"),
        maxRing = cms.int32(7),
        minRing = cms.int32(5),
        rphiRecHits = cms.InputTag("siStripMatchedRecHitsTop","rphiRecHit"),
        useRingSlector = cms.bool(False),
        useSimpleRphiHitsCleaner = cms.bool(True)
    ),
    TIB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHitsTop","matchedRecHit")
    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHitsTop","matchedRecHit")
    ),
    layerList = cms.vstring(
        'MTOB5+MTOB6',
        'MTOB4+MTOB5'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.combinatorialcosmicseedingtripletsP5 = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(

    ),
    FPix = cms.PSet(

    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit")
    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit")
    ),
    TEC = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        maxRing = cms.int32(7),
        minRing = cms.int32(5),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
        useRingSlector = cms.bool(False),
        useSimpleRphiHitsCleaner = cms.bool(True)
    ),
    TIB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit")
    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit")
    ),
    layerList = cms.vstring(
        'MTOB4+MTOB5+MTOB6',
        'MTOB3+MTOB5+MTOB6',
        'MTOB3+MTOB4+MTOB5',
        'TOB2+MTOB4+MTOB5',
        'MTOB3+MTOB4+MTOB6',
        'TOB2+MTOB4+MTOB6'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.combinatorialcosmicseedingtripletsP5Bottom = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(

    ),
    FPix = cms.PSet(

    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        rphiRecHits = cms.InputTag("siStripMatchedRecHitsBottom","rphiRecHit")
    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        rphiRecHits = cms.InputTag("siStripMatchedRecHitsBottom","rphiRecHit")
    ),
    TEC = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHitsBottom","matchedRecHit"),
        maxRing = cms.int32(7),
        minRing = cms.int32(5),
        rphiRecHits = cms.InputTag("siStripMatchedRecHitsBottom","rphiRecHit"),
        useRingSlector = cms.bool(False),
        useSimpleRphiHitsCleaner = cms.bool(True)
    ),
    TIB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHitsBottom","matchedRecHit")
    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHitsBottom","matchedRecHit")
    ),
    layerList = cms.vstring(
        'MTOB4+MTOB5+MTOB6',
        'MTOB3+MTOB5+MTOB6',
        'MTOB3+MTOB4+MTOB5',
        'TOB2+MTOB4+MTOB5',
        'MTOB3+MTOB4+MTOB6',
        'TOB2+MTOB4+MTOB6'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.combinatorialcosmicseedingtripletsP5Top = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(

    ),
    FPix = cms.PSet(

    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        rphiRecHits = cms.InputTag("siStripMatchedRecHitsTop","rphiRecHit")
    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        rphiRecHits = cms.InputTag("siStripMatchedRecHitsTop","rphiRecHit")
    ),
    TEC = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHitsTop","matchedRecHit"),
        maxRing = cms.int32(7),
        minRing = cms.int32(5),
        rphiRecHits = cms.InputTag("siStripMatchedRecHitsTop","rphiRecHit"),
        useRingSlector = cms.bool(False),
        useSimpleRphiHitsCleaner = cms.bool(True)
    ),
    TIB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHitsTop","matchedRecHit")
    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHitsTop","matchedRecHit")
    ),
    layerList = cms.vstring(
        'MTOB4+MTOB5+MTOB6',
        'MTOB3+MTOB5+MTOB6',
        'MTOB3+MTOB4+MTOB5',
        'TOB2+MTOB4+MTOB5',
        'MTOB3+MTOB4+MTOB6',
        'TOB2+MTOB4+MTOB6'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.combinatoricRecoTaus = cms.EDProducer("RecoTauProducer",
    buildNullTaus = cms.bool(False),
    builders = cms.VPSet(cms.PSet(
        decayModes = cms.VPSet(
            cms.PSet(
                maxPiZeros = cms.uint32(0),
                maxTracks = cms.uint32(6),
                nCharged = cms.uint32(1),
                nPiZeros = cms.uint32(0)
            ),
            cms.PSet(
                maxPiZeros = cms.uint32(6),
                maxTracks = cms.uint32(6),
                nCharged = cms.uint32(1),
                nPiZeros = cms.uint32(1)
            ),
            cms.PSet(
                maxPiZeros = cms.uint32(5),
                maxTracks = cms.uint32(6),
                nCharged = cms.uint32(1),
                nPiZeros = cms.uint32(2)
            ),
            cms.PSet(
                maxPiZeros = cms.uint32(0),
                maxTracks = cms.uint32(6),
                nCharged = cms.uint32(2),
                nPiZeros = cms.uint32(0)
            ),
            cms.PSet(
                maxPiZeros = cms.uint32(3),
                maxTracks = cms.uint32(6),
                nCharged = cms.uint32(2),
                nPiZeros = cms.uint32(1)
            ),
            cms.PSet(
                maxPiZeros = cms.uint32(0),
                maxTracks = cms.uint32(6),
                nCharged = cms.uint32(3),
                nPiZeros = cms.uint32(0)
            ),
            cms.PSet(
                maxPiZeros = cms.uint32(3),
                maxTracks = cms.uint32(6),
                nCharged = cms.uint32(3),
                nPiZeros = cms.uint32(1)
            )
        ),
        isolationConeSize = cms.double(0.5),
        minAbsPhotonSumPt_insideSignalCone = cms.double(2.5),
        minAbsPhotonSumPt_outsideSignalCone = cms.double(1000000000.0),
        minRelPhotonSumPt_insideSignalCone = cms.double(0.1),
        minRelPhotonSumPt_outsideSignalCone = cms.double(1000000000.0),
        name = cms.string('combinatoric'),
        pfCandSrc = cms.InputTag("particleFlow"),
        plugin = cms.string('RecoTauBuilderCombinatoricPlugin'),
        qualityCuts = cms.PSet(
            isolationQualityCuts = cms.PSet(
                maxDeltaZ = cms.double(0.2),
                maxDeltaZToLeadTrack = cms.double(-1.0),
                maxTrackChi2 = cms.double(100.0),
                maxTransverseImpactParameter = cms.double(0.03),
                minGammaEt = cms.double(1.5),
                minTrackHits = cms.uint32(8),
                minTrackPixelHits = cms.uint32(0),
                minTrackPt = cms.double(1.0),
                minTrackVertexWeight = cms.double(-1.0)
            ),
            leadingTrkOrPFCandOption = cms.string('leadPFCand'),
            primaryVertexSrc = cms.InputTag("offlinePrimaryVertices"),
            pvFindingAlgo = cms.string('closestInDeltaZ'),
            recoverLeadingTrk = cms.bool(False),
            signalQualityCuts = cms.PSet(
                maxDeltaZ = cms.double(0.4),
                maxDeltaZToLeadTrack = cms.double(-1.0),
                maxTrackChi2 = cms.double(100.0),
                maxTransverseImpactParameter = cms.double(0.1),
                minGammaEt = cms.double(1.0),
                minNeutralHadronEt = cms.double(30.0),
                minTrackHits = cms.uint32(3),
                minTrackPixelHits = cms.uint32(0),
                minTrackPt = cms.double(0.5),
                minTrackVertexWeight = cms.double(-1.0)
            ),
            vertexTrackFiltering = cms.bool(False),
            vxAssocQualityCuts = cms.PSet(
                maxTrackChi2 = cms.double(100.0),
                maxTransverseImpactParameter = cms.double(0.1),
                minGammaEt = cms.double(1.0),
                minTrackHits = cms.uint32(3),
                minTrackPixelHits = cms.uint32(0),
                minTrackPt = cms.double(0.5),
                minTrackVertexWeight = cms.double(-1.0)
            )
        ),
        signalConeSize = cms.string('max(min(0.1, 3.0/pt()), 0.05)'),
        verbosity = cms.int32(0)
    )),
    chargedHadronSrc = cms.InputTag("ak4PFJetsRecoTauChargedHadrons"),
    jetRegionSrc = cms.InputTag("recoTauAK4PFJets08Region"),
    jetSrc = cms.InputTag("ak4PFJets"),
    maxJetAbsEta = cms.double(2.7),
    minJetPt = cms.double(14.0),
    modifiers = cms.VPSet(
        cms.PSet(
            name = cms.string('sipt'),
            plugin = cms.string('RecoTauImpactParameterSignificancePlugin'),
            qualityCuts = cms.PSet(
                isolationQualityCuts = cms.PSet(
                    maxDeltaZ = cms.double(0.2),
                    maxDeltaZToLeadTrack = cms.double(-1.0),
                    maxTrackChi2 = cms.double(100.0),
                    maxTransverseImpactParameter = cms.double(0.03),
                    minGammaEt = cms.double(1.5),
                    minTrackHits = cms.uint32(8),
                    minTrackPixelHits = cms.uint32(0),
                    minTrackPt = cms.double(1.0),
                    minTrackVertexWeight = cms.double(-1.0)
                ),
                leadingTrkOrPFCandOption = cms.string('leadPFCand'),
                primaryVertexSrc = cms.InputTag("offlinePrimaryVertices"),
                pvFindingAlgo = cms.string('closestInDeltaZ'),
                recoverLeadingTrk = cms.bool(False),
                signalQualityCuts = cms.PSet(
                    maxDeltaZ = cms.double(0.4),
                    maxDeltaZToLeadTrack = cms.double(-1.0),
                    maxTrackChi2 = cms.double(100.0),
                    maxTransverseImpactParameter = cms.double(0.1),
                    minGammaEt = cms.double(1.0),
                    minNeutralHadronEt = cms.double(30.0),
                    minTrackHits = cms.uint32(3),
                    minTrackPixelHits = cms.uint32(0),
                    minTrackPt = cms.double(0.5),
                    minTrackVertexWeight = cms.double(-1.0)
                ),
                vertexTrackFiltering = cms.bool(False),
                vxAssocQualityCuts = cms.PSet(
                    maxTrackChi2 = cms.double(100.0),
                    maxTransverseImpactParameter = cms.double(0.1),
                    minGammaEt = cms.double(1.0),
                    minTrackHits = cms.uint32(3),
                    minTrackPixelHits = cms.uint32(0),
                    minTrackPt = cms.double(0.5),
                    minTrackVertexWeight = cms.double(-1.0)
                )
            )
        ),
        cms.PSet(
            DataType = cms.string('AOD'),
            EcalStripSumE_deltaEta = cms.double(0.03),
            EcalStripSumE_deltaPhiOverQ_maxValue = cms.double(0.5),
            EcalStripSumE_deltaPhiOverQ_minValue = cms.double(-0.1),
            EcalStripSumE_minClusEnergy = cms.double(0.1),
            ElecPreIDLeadTkMatch_maxDR = cms.double(0.01),
            ElectronPreIDProducer = cms.InputTag("elecpreid"),
            maximumForElectrionPreIDOutput = cms.double(-0.1),
            name = cms.string('elec_rej'),
            plugin = cms.string('RecoTauElectronRejectionPlugin')
        ),
        cms.PSet(
            dRaddNeutralHadron = cms.double(0.12),
            dRaddPhoton = cms.double(-1.0),
            minGammaEt = cms.double(10.0),
            minNeutralHadronEt = cms.double(50.0),
            name = cms.string('tau_en_reconstruction'),
            plugin = cms.string('PFRecoTauEnergyAlgorithmPlugin'),
            verbosity = cms.int32(0)
        ),
        cms.PSet(
            name = cms.string('tau_lost_tracks'),
            plugin = cms.string('PFRecoTauLostTrackPlugin'),
            trackSrc = cms.InputTag("generalTracks"),
            verbosity = cms.int32(0)
        ),
        cms.PSet(
            name = cms.string('tau_mass'),
            plugin = cms.string('PFRecoTauMassPlugin'),
            verbosity = cms.int32(0)
        ),
        cms.PSet(
            name = cms.string('TTIworkaround'),
            pfTauTagInfoSrc = cms.InputTag("pfRecoTauTagInfoProducer"),
            plugin = cms.string('RecoTauTagInfoWorkaroundModifer')
        )
    ),
    outputSelection = cms.string('leadChargedHadrCand().isNonnull()'),
    piZeroSrc = cms.InputTag("ak4PFJetsLegacyHPSPiZeros")
)


process.combinedP5SeedsForCTF = cms.EDProducer("SeedCombiner",
    PairCollection = cms.InputTag("combinatorialcosmicseedfinderP5"),
    TripletCollection = cms.InputTag("simpleCosmicBONSeeds"),
    seedCollections = cms.VInputTag("combinatorialcosmicseedfinderP5", "simpleCosmicBONSeeds")
)


process.combinedP5SeedsForCTFBottom = cms.EDProducer("SeedCombiner",
    seedCollections = cms.VInputTag("combinatorialcosmicseedfinderP5Bottom", "simpleCosmicBONSeedsBottom")
)


process.combinedP5SeedsForCTFTop = cms.EDProducer("SeedCombiner",
    seedCollections = cms.VInputTag("combinatorialcosmicseedfinderP5Top", "simpleCosmicBONSeedsTop")
)


process.conv2Clusters = cms.EDProducer("TrackClusterRemover",
    TrackQuality = cms.string('highPurity'),
    maxChi2 = cms.double(30),
    mightGet = cms.optional.untracked.vstring,
    minNumberOfLayersWithMeasBeforeFiltering = cms.int32(0),
    oldClusterRemovalInfo = cms.InputTag("convClusters"),
    overrideTrkQuals = cms.InputTag("convStepSelector","convStep"),
    pixelClusters = cms.InputTag("siPixelClusters"),
    stripClusters = cms.InputTag("siStripClusters"),
    trackClassifier = cms.InputTag("","QualityMasks"),
    trajectories = cms.InputTag("convStepTracks")
)


process.conv2LayerPairs = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(
        HitProducer = cms.string('siPixelRecHits'),
        TTRHBuilder = cms.string('WithTrackAngle'),
        skipClusters = cms.InputTag("conv2Clusters")
    ),
    FPix = cms.PSet(
        HitProducer = cms.string('siPixelRecHits'),
        TTRHBuilder = cms.string('WithTrackAngle'),
        skipClusters = cms.InputTag("conv2Clusters")
    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTight')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        maxRing = cms.int32(7),
        minRing = cms.int32(1),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHitUnmatched"),
        skipClusters = cms.InputTag("conv2Clusters"),
        stereoRecHits = cms.InputTag("siStripMatchedRecHits","stereoRecHitUnmatched"),
        useRingSlector = cms.bool(True),
        useSimpleRphiHitsCleaner = cms.bool(False)
    ),
    TIB = cms.PSet(

    ),
    TIB1 = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTight')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        skipClusters = cms.InputTag("conv2Clusters")
    ),
    TIB2 = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTight')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        skipClusters = cms.InputTag("conv2Clusters")
    ),
    TIB3 = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTight')
        ),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
        skipClusters = cms.InputTag("conv2Clusters")
    ),
    TIB4 = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTight')
        ),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
        skipClusters = cms.InputTag("conv2Clusters")
    ),
    TID = cms.PSet(

    ),
    TID1 = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTight')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        maxRing = cms.int32(2),
        minRing = cms.int32(1),
        skipClusters = cms.InputTag("conv2Clusters"),
        useRingSlector = cms.bool(True),
        useSimpleRphiHitsCleaner = cms.bool(False)
    ),
    TID2 = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTight')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        maxRing = cms.int32(2),
        minRing = cms.int32(1),
        skipClusters = cms.InputTag("conv2Clusters"),
        useRingSlector = cms.bool(True),
        useSimpleRphiHitsCleaner = cms.bool(False)
    ),
    TID3 = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTight')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        maxRing = cms.int32(2),
        minRing = cms.int32(1),
        skipClusters = cms.InputTag("conv2Clusters"),
        useRingSlector = cms.bool(True),
        useSimpleRphiHitsCleaner = cms.bool(False)
    ),
    TOB = cms.PSet(

    ),
    TOB1 = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTight')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        skipClusters = cms.InputTag("conv2Clusters")
    ),
    TOB2 = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTight')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        skipClusters = cms.InputTag("conv2Clusters")
    ),
    TOB3 = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTight')
        ),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
        skipClusters = cms.InputTag("conv2Clusters")
    ),
    TOB4 = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTight')
        ),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
        skipClusters = cms.InputTag("conv2Clusters")
    ),
    TOB5 = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTight')
        ),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
        skipClusters = cms.InputTag("conv2Clusters")
    ),
    TOB6 = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTight')
        ),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
        skipClusters = cms.InputTag("conv2Clusters")
    ),
    layerList = cms.vstring(
        'BPix1+BPix2',
        'BPix2+BPix3',
        'BPix2+FPix1_pos',
        'BPix2+FPix1_neg',
        'BPix2+FPix2_pos',
        'BPix2+FPix2_neg',
        'FPix1_pos+FPix2_pos',
        'FPix1_neg+FPix2_neg',
        'BPix3+TIB1',
        'TIB1+TID1_pos',
        'TIB1+TID1_neg',
        'TIB1+TID2_pos',
        'TIB1+TID2_neg',
        'TIB1+TIB2',
        'TIB2+TID1_pos',
        'TIB2+TID1_neg',
        'TIB2+TID2_pos',
        'TIB2+TID2_neg',
        'TIB2+TIB3',
        'TIB3+TIB4',
        'TIB3+TID1_pos',
        'TIB3+TID1_neg',
        'TIB4+TOB1',
        'TOB1+TOB2',
        'TOB1+TEC1_pos',
        'TOB1+TEC1_neg',
        'TOB2+TOB3',
        'TOB2+TEC1_pos',
        'TOB2+TEC1_neg',
        'TOB3+TOB4',
        'TOB3+TEC1_pos',
        'TOB3+TEC1_neg',
        'TOB4+TOB5',
        'TOB5+TOB6',
        'TID1_pos+TID2_pos',
        'TID2_pos+TID3_pos',
        'TID3_pos+TEC1_pos',
        'TID1_neg+TID2_neg',
        'TID2_neg+TID3_neg',
        'TID3_neg+TEC1_neg',
        'TEC1_pos+TEC2_pos',
        'TEC2_pos+TEC3_pos',
        'TEC3_pos+TEC4_pos',
        'TEC4_pos+TEC5_pos',
        'TEC5_pos+TEC6_pos',
        'TEC6_pos+TEC7_pos',
        'TEC7_pos+TEC8_pos',
        'TEC1_neg+TEC2_neg',
        'TEC2_neg+TEC3_neg',
        'TEC3_neg+TEC4_neg',
        'TEC4_neg+TEC5_neg',
        'TEC5_neg+TEC6_neg',
        'TEC6_neg+TEC7_neg',
        'TEC7_neg+TEC8_neg'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.conv2StepSelector = cms.EDProducer("MultiTrackSelector",
    beamspot = cms.InputTag("offlineBeamSpot"),
    src = cms.InputTag("conv2StepTracks"),
    trackSelectors = cms.VPSet(
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(False),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(3.0),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(5.0, 8.0),
            d0_par2 = cms.vdouble(5.0, 8.0),
            dz_par1 = cms.vdouble(5.0, 8.0),
            dz_par2 = cms.vdouble(5.0, 8.0),
            keepAllTracks = cms.bool(False),
            maxNumberLostLayers = cms.uint32(1),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(1),
            minNumberLayers = cms.uint32(3),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(0),
            nSigmaZ = cms.double(4.0),
            name = cms.string('conv2StepLoose'),
            preFilterName = cms.string(''),
            qualityBit = cms.string('loose'),
            res_par = cms.vdouble(0.003, 0.001),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        ),
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(True),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(2.5),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(5.0, 8.0),
            d0_par2 = cms.vdouble(5.0, 8.0),
            dz_par1 = cms.vdouble(5.0, 8.0),
            dz_par2 = cms.vdouble(5.0, 8.0),
            keepAllTracks = cms.bool(True),
            maxNumberLostLayers = cms.uint32(1),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(1),
            minNumberLayers = cms.uint32(3),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(0),
            nSigmaZ = cms.double(4.0),
            name = cms.string('conv2StepTight'),
            preFilterName = cms.string('conv2StepLoose'),
            qualityBit = cms.string('tight'),
            res_par = cms.vdouble(0.003, 0.001),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        ),
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(True),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(2.0),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(5.0, 8.0),
            d0_par2 = cms.vdouble(5.0, 8.0),
            dz_par1 = cms.vdouble(5.0, 8.0),
            dz_par2 = cms.vdouble(5.0, 8.0),
            keepAllTracks = cms.bool(True),
            maxNumberLostLayers = cms.uint32(1),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(1),
            minNumberLayers = cms.uint32(3),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(0),
            nSigmaZ = cms.double(4.0),
            name = cms.string('conv2Step'),
            preFilterName = cms.string('conv2StepTight'),
            qualityBit = cms.string('highPurity'),
            res_par = cms.vdouble(0.003, 0.001),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        )
    ),
    useVertices = cms.bool(True),
    useVtxError = cms.bool(False),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.conv2StepTracks = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('conversionStep'),
    Fitter = cms.string('conv2StepFitterSmoother'),
    GeometricInnerState = cms.bool(False),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    Propagator = cms.string('RungeKuttaTrackerPropagator'),
    SimpleMagneticField = cms.string(''),
    TTRHBuilder = cms.string('WithAngleAndTemplate'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("conv2TrackCandidates"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(False)
)


process.conv2TrackCandidates = cms.EDProducer("CkfTrackCandidateMaker",
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    RedundantSeedCleaner = cms.string('CachingSeedCleanerBySharedInput'),
    TrajectoryBuilderPSet = cms.PSet(
        refToPSet_ = cms.string('conv2CkfTrajectoryBuilder')
    ),
    TrajectoryCleaner = cms.string('TrajectoryCleanerBySharedHits'),
    TransientInitialStateEstimatorParameters = cms.PSet(
        numberMeasurementsForFit = cms.int32(4),
        propagatorAlongTISE = cms.string('PropagatorWithMaterial'),
        propagatorOppositeTISE = cms.string('PropagatorWithMaterialOpposite')
    ),
    cleanTrajectoryAfterInOut = cms.bool(True),
    clustersToSkip = cms.InputTag(""),
    doSeedingRegionRebuilding = cms.bool(True),
    maxNSeeds = cms.uint32(500000),
    maxSeedsBeforeCleaning = cms.uint32(5000),
    mightGet = cms.optional.untracked.vstring,
    numHitsForSeedCleaner = cms.int32(4),
    onlyPixelHitsForSeedCleaner = cms.bool(False),
    phase2clustersToSkip = cms.InputTag(""),
    reverseTrajectories = cms.bool(False),
    src = cms.InputTag("photonConvTrajSeedFromQuadruplets","conv2SeedCandidates"),
    useHitsSplitting = cms.bool(True)
)


process.convClusters = cms.EDProducer("TrackClusterRemover",
    TrackQuality = cms.string('highPurity'),
    maxChi2 = cms.double(30),
    mightGet = cms.optional.untracked.vstring,
    minNumberOfLayersWithMeasBeforeFiltering = cms.int32(0),
    oldClusterRemovalInfo = cms.InputTag("tobTecStepClusters"),
    overrideTrkQuals = cms.InputTag(""),
    pixelClusters = cms.InputTag("siPixelClusters"),
    stripClusters = cms.InputTag("siStripClusters"),
    trackClassifier = cms.InputTag("tobTecStep","QualityMasks"),
    trajectories = cms.InputTag("tobTecStepTracks")
)


process.convLayerPairs = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(
        HitProducer = cms.string('siPixelRecHits'),
        TTRHBuilder = cms.string('WithTrackAngle'),
        skipClusters = cms.InputTag("convClusters")
    ),
    FPix = cms.PSet(
        HitProducer = cms.string('siPixelRecHits'),
        TTRHBuilder = cms.string('WithTrackAngle'),
        skipClusters = cms.InputTag("convClusters")
    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTight')
        ),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
        skipClusters = cms.InputTag("convClusters")
    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTight')
        ),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
        skipClusters = cms.InputTag("convClusters")
    ),
    TEC = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTight')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        maxRing = cms.int32(7),
        minRing = cms.int32(1),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHitUnmatched"),
        skipClusters = cms.InputTag("convClusters"),
        stereoRecHits = cms.InputTag("siStripMatchedRecHits","stereoRecHitUnmatched"),
        useRingSlector = cms.bool(True),
        useSimpleRphiHitsCleaner = cms.bool(False)
    ),
    TIB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTight')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        skipClusters = cms.InputTag("convClusters")
    ),
    TID = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTight')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        maxRing = cms.int32(2),
        minRing = cms.int32(1),
        skipClusters = cms.InputTag("convClusters"),
        useRingSlector = cms.bool(True),
        useSimpleRphiHitsCleaner = cms.bool(False)
    ),
    TOB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTight')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        skipClusters = cms.InputTag("convClusters")
    ),
    layerList = cms.vstring(
        'BPix1+BPix2',
        'BPix2+BPix3',
        'BPix3+BPix4',
        'BPix2+FPix1_pos',
        'BPix2+FPix1_neg',
        'BPix2+FPix2_pos',
        'BPix2+FPix2_neg',
        'BPix3+FPix1_pos',
        'BPix3+FPix1_neg',
        'FPix1_pos+FPix2_pos',
        'FPix1_neg+FPix2_neg',
        'FPix2_pos+FPix3_pos',
        'FPix2_neg+FPix3_neg',
        'BPix3+TIB1',
        'BPix4+TIB1',
        'BPix4+TIB2',
        'TIB1+TID1_pos',
        'TIB1+TID1_neg',
        'TIB1+TID2_pos',
        'TIB1+TID2_neg',
        'TIB1+TIB2',
        'TIB1+MTIB3',
        'TIB2+TID1_pos',
        'TIB2+TID1_neg',
        'TIB2+TID2_pos',
        'TIB2+TID2_neg',
        'TIB2+MTIB3',
        'TIB2+MTIB4',
        'MTIB3+MTIB4',
        'MTIB3+TOB1',
        'MTIB3+TID1_pos',
        'MTIB3+TID1_neg',
        'MTIB4+TOB1',
        'MTIB4+TOB2',
        'TOB1+TOB2',
        'TOB1+MTOB3',
        'TOB1+TEC1_pos',
        'TOB1+TEC1_neg',
        'TOB2+MTOB3',
        'TOB2+MTOB4',
        'TOB2+TEC1_pos',
        'TOB2+TEC1_neg',
        'TID1_pos+TID2_pos',
        'TID2_pos+TID3_pos',
        'TID3_pos+TEC1_pos',
        'TID1_neg+TID2_neg',
        'TID2_neg+TID3_neg',
        'TID3_neg+TEC1_neg',
        'TEC1_pos+TEC2_pos',
        'TEC2_pos+TEC3_pos',
        'TEC3_pos+TEC4_pos',
        'TEC4_pos+TEC5_pos',
        'TEC5_pos+TEC6_pos',
        'TEC6_pos+TEC7_pos',
        'TEC7_pos+TEC8_pos',
        'TEC1_neg+TEC2_neg',
        'TEC2_neg+TEC3_neg',
        'TEC3_neg+TEC4_neg',
        'TEC4_neg+TEC5_neg',
        'TEC5_neg+TEC6_neg',
        'TEC6_neg+TEC7_neg',
        'TEC7_neg+TEC8_neg'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.convStepSelector = cms.EDProducer("MultiTrackSelector",
    beamspot = cms.InputTag("offlineBeamSpot"),
    src = cms.InputTag("convStepTracks"),
    trackSelectors = cms.VPSet(
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(False),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(3.0),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(5.0, 8.0),
            d0_par2 = cms.vdouble(5.0, 8.0),
            dz_par1 = cms.vdouble(5.0, 8.0),
            dz_par2 = cms.vdouble(5.0, 8.0),
            keepAllTracks = cms.bool(False),
            maxNumberLostLayers = cms.uint32(1),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(1),
            minNumberLayers = cms.uint32(3),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(0),
            nSigmaZ = cms.double(4.0),
            name = cms.string('convStepLoose'),
            preFilterName = cms.string(''),
            qualityBit = cms.string('loose'),
            res_par = cms.vdouble(0.003, 0.001),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        ),
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(True),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(2.5),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(5.0, 8.0),
            d0_par2 = cms.vdouble(5.0, 8.0),
            dz_par1 = cms.vdouble(5.0, 8.0),
            dz_par2 = cms.vdouble(5.0, 8.0),
            keepAllTracks = cms.bool(True),
            maxNumberLostLayers = cms.uint32(1),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(1),
            minNumberLayers = cms.uint32(3),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(0),
            nSigmaZ = cms.double(4.0),
            name = cms.string('convStepTight'),
            preFilterName = cms.string('convStepLoose'),
            qualityBit = cms.string('tight'),
            res_par = cms.vdouble(0.003, 0.001),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        ),
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(True),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(2.0),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(5.0, 8.0),
            d0_par2 = cms.vdouble(5.0, 8.0),
            dz_par1 = cms.vdouble(5.0, 8.0),
            dz_par2 = cms.vdouble(5.0, 8.0),
            keepAllTracks = cms.bool(True),
            maxNumberLostLayers = cms.uint32(1),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(1),
            minNumberLayers = cms.uint32(3),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(0),
            nSigmaZ = cms.double(4.0),
            name = cms.string('convStep'),
            preFilterName = cms.string('convStepTight'),
            qualityBit = cms.string('highPurity'),
            res_par = cms.vdouble(0.003, 0.001),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        )
    ),
    useVertices = cms.bool(True),
    useVtxError = cms.bool(False),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.convStepTracks = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('conversionStep'),
    Fitter = cms.string('convStepFitterSmoother'),
    GeometricInnerState = cms.bool(False),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    Propagator = cms.string('RungeKuttaTrackerPropagator'),
    SimpleMagneticField = cms.string(''),
    TTRHBuilder = cms.string('WithAngleAndTemplate'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("convTrackCandidates"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(False)
)


process.convTrackCandidates = cms.EDProducer("CkfTrackCandidateMaker",
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    RedundantSeedCleaner = cms.string('CachingSeedCleanerBySharedInput'),
    TrajectoryBuilderPSet = cms.PSet(
        refToPSet_ = cms.string('convCkfTrajectoryBuilder')
    ),
    TrajectoryCleaner = cms.string('TrajectoryCleanerBySharedHits'),
    TransientInitialStateEstimatorParameters = cms.PSet(
        numberMeasurementsForFit = cms.int32(4),
        propagatorAlongTISE = cms.string('PropagatorWithMaterial'),
        propagatorOppositeTISE = cms.string('PropagatorWithMaterialOpposite')
    ),
    cleanTrajectoryAfterInOut = cms.bool(True),
    clustersToSkip = cms.InputTag("convClusters"),
    doSeedingRegionRebuilding = cms.bool(True),
    maxNSeeds = cms.uint32(500000),
    maxSeedsBeforeCleaning = cms.uint32(5000),
    mightGet = cms.optional.untracked.vstring,
    numHitsForSeedCleaner = cms.int32(4),
    onlyPixelHitsForSeedCleaner = cms.bool(False),
    phase2clustersToSkip = cms.InputTag(""),
    reverseTrajectories = cms.bool(False),
    src = cms.InputTag("photonConvTrajSeedFromSingleLeg","convSeedCandidates"),
    useHitsSplitting = cms.bool(True)
)


process.conversionStepTracks = cms.EDProducer("TrackListMerger",
    Epsilon = cms.double(-0.001),
    FoundHitBonus = cms.double(5.0),
    LostHitPenalty = cms.double(5.0),
    MaxNormalizedChisq = cms.double(1000.0),
    MinFound = cms.int32(3),
    MinPT = cms.double(0.05),
    ShareFrac = cms.double(0.19),
    TrackProducers = cms.VInputTag("convStepTracks"),
    allowFirstHitShare = cms.bool(True),
    copyExtras = cms.untracked.bool(True),
    copyMVA = cms.bool(True),
    hasSelector = cms.vint32(1),
    indivShareFrac = cms.vdouble(1.0, 1.0),
    makeReKeyedSeeds = cms.untracked.bool(False),
    newQuality = cms.string('confirmed'),
    selectedTrackQuals = cms.VInputTag("convStepSelector:convStep"),
    setsToMerge = cms.VPSet(cms.PSet(
        pQual = cms.bool(True),
        tLists = cms.vint32(1)
    )),
    trackAlgoPriorityOrder = cms.string('trackAlgoPriorityOrder'),
    writeOnlyTrkQuals = cms.bool(False)
)


process.cosmicCandidateFinder = cms.EDProducer("CosmicTrackFinder",
    Chi2Cut = cms.double(30.0),
    GeometricStructure = cms.untracked.string('MTCC'),
    HitProducer = cms.string('siStripRecHits'),
    MinHits = cms.int32(4),
    TTRHBuilder = cms.string('WithTrackAngle'),
    cosmicSeeds = cms.InputTag("cosmicseedfinder"),
    debug = cms.untracked.bool(True),
    matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
    pixelRecHits = cms.InputTag("siPixelRecHits"),
    rphirecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
    stereorecHits = cms.InputTag("siStripMatchedRecHits","stereoRecHit"),
    useHitsSplitting = cms.bool(True)
)


process.cosmicCandidateFinderP5 = cms.EDProducer("CosmicTrackFinder",
    Chi2Cut = cms.double(30.0),
    GeometricStructure = cms.untracked.string('STANDARD'),
    HitProducer = cms.string('siStripRecHits'),
    MinHits = cms.int32(4),
    TTRHBuilder = cms.string('WithTrackAngle'),
    cosmicSeeds = cms.InputTag("cosmicseedfinderP5"),
    debug = cms.untracked.bool(True),
    matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
    pixelRecHits = cms.InputTag("siPixelRecHits"),
    rphirecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
    stereorecHits = cms.InputTag("siStripMatchedRecHits","stereoRecHit"),
    useHitsSplitting = cms.bool(True)
)


process.cosmicCandidateFinderP5Bottom = cms.EDProducer("CosmicTrackFinder",
    Chi2Cut = cms.double(30.0),
    GeometricStructure = cms.untracked.string('STANDARD'),
    HitProducer = cms.string('siStripRecHitsBottom'),
    MinHits = cms.int32(4),
    TTRHBuilder = cms.string('WithTrackAngle'),
    cosmicSeeds = cms.InputTag("cosmicseedfinderP5Bottom"),
    debug = cms.untracked.bool(True),
    matchedRecHits = cms.InputTag("siStripMatchedRecHitsBottom","matchedRecHit"),
    pixelRecHits = cms.InputTag("siPixelRecHitsBottom"),
    rphirecHits = cms.InputTag("siStripMatchedRecHitsBottom","rphiRecHit"),
    stereorecHits = cms.InputTag("siStripMatchedRecHitsBottom","stereoRecHit"),
    useHitsSplitting = cms.bool(True)
)


process.cosmicCandidateFinderP5Top = cms.EDProducer("CosmicTrackFinder",
    Chi2Cut = cms.double(30.0),
    GeometricStructure = cms.untracked.string('STANDARD'),
    HitProducer = cms.string('siStripRecHitsTop'),
    MinHits = cms.int32(4),
    TTRHBuilder = cms.string('WithTrackAngle'),
    cosmicSeeds = cms.InputTag("cosmicseedfinderP5Top"),
    debug = cms.untracked.bool(True),
    matchedRecHits = cms.InputTag("siStripMatchedRecHitsTop","matchedRecHit"),
    pixelRecHits = cms.InputTag("siPixelRecHitsTop"),
    rphirecHits = cms.InputTag("siStripMatchedRecHitsTop","rphiRecHit"),
    stereorecHits = cms.InputTag("siStripMatchedRecHitsTop","stereoRecHit"),
    useHitsSplitting = cms.bool(True)
)


process.cosmicDCCkfTrackCandidates = cms.EDProducer("CkfTrackCandidateMaker",
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('CosmicNavigationSchool'),
    RedundantSeedCleaner = cms.string('CachingSeedCleanerBySharedInput'),
    TrajectoryBuilderPSet = cms.PSet(
        refToPSet_ = cms.string('GroupedCkfTrajectoryBuilderCDC')
    ),
    TrajectoryCleaner = cms.string('TrajectoryCleanerBySharedHits'),
    TransientInitialStateEstimatorParameters = cms.PSet(
        numberMeasurementsForFit = cms.int32(4),
        propagatorAlongTISE = cms.string('PropagatorWithMaterial'),
        propagatorOppositeTISE = cms.string('PropagatorWithMaterialOpposite')
    ),
    cleanTrajectoryAfterInOut = cms.bool(True),
    clustersToSkip = cms.InputTag(""),
    doSeedingRegionRebuilding = cms.bool(True),
    maxNSeeds = cms.uint32(500000),
    maxSeedsBeforeCleaning = cms.uint32(5000),
    mightGet = cms.optional.untracked.vstring,
    numHitsForSeedCleaner = cms.int32(4),
    onlyPixelHitsForSeedCleaner = cms.bool(False),
    phase2clustersToSkip = cms.InputTag(""),
    reverseTrajectories = cms.bool(False),
    src = cms.InputTag("cosmicDCSeeds"),
    useHitsSplitting = cms.bool(True)
)


process.cosmicDCSeeds = cms.EDProducer("OutsideInMuonSeeder",
    cut = cms.string('p > 3 && abs(eta)<1.6 && phi<0'),
    debug = cms.untracked.bool(False),
    errorRescaleFactor = cms.double(2.0),
    fromVertex = cms.bool(False),
    hitCollector = cms.string('hitCollectorForCosmicDCSeeds'),
    hitsToTry = cms.int32(3),
    layersToTry = cms.int32(3),
    maxEtaForTOB = cms.double(2.5),
    minEtaForTEC = cms.double(0.0),
    muonPropagator = cms.string('SteppingHelixPropagatorAlong'),
    src = cms.InputTag("muonsFromCosmics"),
    trackerPropagator = cms.string('PropagatorWithMaterial')
)


process.cosmicDCTracks = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('ctf'),
    Fitter = cms.string('FittingSmootherRKP5'),
    GeometricInnerState = cms.bool(True),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    Propagator = cms.string('RungeKuttaTrackerPropagator'),
    SimpleMagneticField = cms.string(''),
    TTRHBuilder = cms.string('WithAngleAndTemplate'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("cosmicDCCkfTrackCandidates"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(False)
)


process.cosmicTrackSplitter = cms.EDProducer("CosmicTrackSplitter",
    detsToIgnore = cms.vuint32(),
    dxyCut = cms.double(9999.0),
    dzCut = cms.double(9999.0),
    excludePixelHits = cms.bool(False),
    minimumHits = cms.uint32(6),
    replaceWithInactiveHits = cms.bool(False),
    stripAllInvalidHits = cms.bool(False),
    stripBackInvalidHits = cms.bool(True),
    stripFrontInvalidHits = cms.bool(True),
    tjTkAssociationMapTag = cms.InputTag("cosmictrackfinderP5"),
    tracks = cms.InputTag("cosmictrackfinderP5")
)


process.cosmicTrackSplitting = cms.EDProducer("CosmicTrackSplitter",
    detsToIgnore = cms.vuint32(),
    dxyCut = cms.double(9999.0),
    dzCut = cms.double(9999.0),
    excludePixelHits = cms.bool(False),
    minimumHits = cms.uint32(6),
    replaceWithInactiveHits = cms.bool(False),
    stripAllInvalidHits = cms.bool(False),
    stripBackInvalidHits = cms.bool(True),
    stripFrontInvalidHits = cms.bool(True),
    tjTkAssociationMapTag = cms.InputTag("cosmictrackfinderCosmics"),
    tracks = cms.InputTag("cosmictrackfinderCosmics")
)


process.cosmicseedfinder = cms.EDProducer("CosmicSeedGenerator",
    ClusterCollectionLabel = cms.InputTag("siStripClusters"),
    DontCountDetsAboveNClusters = cms.uint32(20),
    GeometricStructure = cms.untracked.string('STANDARD'),
    HitsForSeeds = cms.untracked.string('pairs'),
    MaxNumberOfCosmicClusters = cms.uint32(300),
    MaxNumberOfPixelClusters = cms.uint32(300),
    NegativeYOnly = cms.bool(False),
    PixelClusterCollectionLabel = cms.InputTag("siPixelClusters"),
    PositiveYOnly = cms.bool(False),
    SeedPt = cms.double(5.0),
    TTRHBuilder = cms.string('WithTrackAngle'),
    doClusterCheck = cms.bool(True),
    matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
    maxSeeds = cms.int32(10000),
    originHalfLength = cms.double(90.0),
    originRadius = cms.double(150.0),
    originZPosition = cms.double(0.0),
    ptMin = cms.double(0.9),
    rphirecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
    stereorecHits = cms.InputTag("siStripMatchedRecHits","stereoRecHit")
)


process.cosmicseedfinderP5 = cms.EDProducer("CosmicSeedGenerator",
    ClusterCollectionLabel = cms.InputTag("siStripClusters"),
    DontCountDetsAboveNClusters = cms.uint32(20),
    GeometricStructure = cms.untracked.string('STANDARD'),
    HitsForSeeds = cms.untracked.string('pairs'),
    MaxNumberOfCosmicClusters = cms.uint32(300),
    MaxNumberOfPixelClusters = cms.uint32(300),
    NegativeYOnly = cms.bool(False),
    PixelClusterCollectionLabel = cms.InputTag("siPixelClusters"),
    PositiveYOnly = cms.bool(False),
    SeedPt = cms.double(5.0),
    TTRHBuilder = cms.string('WithTrackAngle'),
    doClusterCheck = cms.bool(True),
    matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
    maxSeeds = cms.int32(10000),
    originHalfLength = cms.double(90.0),
    originRadius = cms.double(150.0),
    originZPosition = cms.double(0.0),
    ptMin = cms.double(0.9),
    rphirecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
    stereorecHits = cms.InputTag("siStripMatchedRecHits","stereoRecHit")
)


process.cosmicseedfinderP5Bottom = cms.EDProducer("CosmicSeedGenerator",
    ClusterCollectionLabel = cms.InputTag("siStripClustersBottom"),
    DontCountDetsAboveNClusters = cms.uint32(20),
    GeometricStructure = cms.untracked.string('STANDARD'),
    HitsForSeeds = cms.untracked.string('pairs'),
    MaxNumberOfCosmicClusters = cms.uint32(150),
    MaxNumberOfPixelClusters = cms.uint32(300),
    NegativeYOnly = cms.bool(True),
    PixelClusterCollectionLabel = cms.InputTag("siPixelClusters"),
    PositiveYOnly = cms.bool(False),
    SeedPt = cms.double(5.0),
    TTRHBuilder = cms.string('WithTrackAngle'),
    doClusterCheck = cms.bool(True),
    matchedRecHits = cms.InputTag("siStripMatchedRecHitsBottom","matchedRecHit"),
    maxSeeds = cms.int32(10000),
    originHalfLength = cms.double(90.0),
    originRadius = cms.double(150.0),
    originZPosition = cms.double(0.0),
    ptMin = cms.double(0.9),
    rphirecHits = cms.InputTag("siStripMatchedRecHitsBottom","rphiRecHit"),
    stereorecHits = cms.InputTag("siStripMatchedRecHitsBottom","stereoRecHit")
)


process.cosmicseedfinderP5Top = cms.EDProducer("CosmicSeedGenerator",
    ClusterCollectionLabel = cms.InputTag("siStripClustersTop"),
    DontCountDetsAboveNClusters = cms.uint32(20),
    GeometricStructure = cms.untracked.string('STANDARD'),
    HitsForSeeds = cms.untracked.string('pairs'),
    MaxNumberOfCosmicClusters = cms.uint32(150),
    MaxNumberOfPixelClusters = cms.uint32(300),
    NegativeYOnly = cms.bool(False),
    PixelClusterCollectionLabel = cms.InputTag("siPixelClusters"),
    PositiveYOnly = cms.bool(True),
    SeedPt = cms.double(5.0),
    TTRHBuilder = cms.string('WithTrackAngle'),
    doClusterCheck = cms.bool(True),
    matchedRecHits = cms.InputTag("siStripMatchedRecHitsTop","matchedRecHit"),
    maxSeeds = cms.int32(10000),
    originHalfLength = cms.double(90.0),
    originRadius = cms.double(150.0),
    originZPosition = cms.double(0.0),
    ptMin = cms.double(0.9),
    rphirecHits = cms.InputTag("siStripMatchedRecHitsTop","rphiRecHit"),
    stereorecHits = cms.InputTag("siStripMatchedRecHitsTop","stereoRecHit")
)


process.cosmictrackfinderCosmics = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('cosmic'),
    Fitter = cms.string('RKFittingSmoother'),
    GeometricInnerState = cms.bool(False),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    Propagator = cms.string('RungeKuttaTrackerPropagator'),
    SimpleMagneticField = cms.string(''),
    TTRHBuilder = cms.string('WithTrackAngle'),
    TrajectoryInEvent = cms.bool(True),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("cosmicCandidateFinderP5"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(False)
)


process.cosmictrackfinderP5 = cms.EDProducer("CosmicTrackSelector",
    beamspot = cms.InputTag("offlineBeamSpot"),
    chi2n_par = cms.double(10.0),
    copyExtras = cms.untracked.bool(True),
    copyTrajectories = cms.untracked.bool(False),
    keepAllTracks = cms.bool(False),
    maxNumberLostLayers = cms.uint32(999),
    max_d0 = cms.double(110.0),
    max_eta = cms.double(2.0),
    max_z0 = cms.double(300.0),
    minNumber3DLayers = cms.uint32(0),
    minNumberLayers = cms.uint32(0),
    min_nHit = cms.uint32(5),
    min_nPixelHit = cms.uint32(0),
    min_pt = cms.double(1.0),
    qualityBit = cms.string(''),
    src = cms.InputTag("cosmictrackfinderCosmics")
)


process.cosmictrackfinderP5Bottom = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('cosmic'),
    Fitter = cms.string('RKFittingSmoother'),
    GeometricInnerState = cms.bool(False),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    Propagator = cms.string('RungeKuttaTrackerPropagator'),
    SimpleMagneticField = cms.string(''),
    TTRHBuilder = cms.string('WithTrackAngle'),
    TrajectoryInEvent = cms.bool(True),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag("topBottomClusterInfoProducerBottom"),
    src = cms.InputTag("cosmicCandidateFinderP5Bottom"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(False)
)


process.cosmictrackfinderP5Top = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('cosmic'),
    Fitter = cms.string('RKFittingSmoother'),
    GeometricInnerState = cms.bool(False),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    Propagator = cms.string('RungeKuttaTrackerPropagator'),
    SimpleMagneticField = cms.string(''),
    TTRHBuilder = cms.string('WithTrackAngle'),
    TrajectoryInEvent = cms.bool(True),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag("topBottomClusterInfoProducerTop"),
    src = cms.InputTag("cosmicCandidateFinderP5Top"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(False)
)


process.ctfCombinedSeeds = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('undefAlgorithm'),
    Fitter = cms.string('KFFittingSmootherWithOutliersRejectionAndRK'),
    GeometricInnerState = cms.bool(False),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    Propagator = cms.string('RungeKuttaTrackerPropagator'),
    SimpleMagneticField = cms.string(''),
    TTRHBuilder = cms.string('WithAngleAndTemplate'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("ckfTrackCandidatesCombinedSeeds"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(False)
)


process.ctfNoOverlaps = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('undefAlgorithm'),
    Fitter = cms.string('KFFittingSmootherWithOutliersRejectionAndRK'),
    GeometricInnerState = cms.bool(False),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    Propagator = cms.string('RungeKuttaTrackerPropagator'),
    SimpleMagneticField = cms.string(''),
    TTRHBuilder = cms.string('WithAngleAndTemplate'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("ckfTrackCandidatesNoOverlaps"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(False)
)


process.ctfPixelLess = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('undefAlgorithm'),
    Fitter = cms.string('RKFittingSmoother'),
    GeometricInnerState = cms.bool(False),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    Propagator = cms.string('RungeKuttaTrackerPropagator'),
    SimpleMagneticField = cms.string(''),
    TTRHBuilder = cms.string('WithAngleAndTemplate'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("ckfTrackCandidatesPixelLess"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(False)
)


process.ctfWithMaterialTracks = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('undefAlgorithm'),
    Fitter = cms.string('KFFittingSmootherWithOutliersRejectionAndRK'),
    GeometricInnerState = cms.bool(False),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    Propagator = cms.string('RungeKuttaTrackerPropagator'),
    SimpleMagneticField = cms.string(''),
    TTRHBuilder = cms.string('WithAngleAndTemplate'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("ckfTrackCandidates"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(False)
)


process.ctfWithMaterialTracksCosmics = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('ctf'),
    Fitter = cms.string('FittingSmootherRKP5'),
    GeometricInnerState = cms.bool(True),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    Propagator = cms.string('RungeKuttaTrackerPropagator'),
    SimpleMagneticField = cms.string(''),
    TTRHBuilder = cms.string('WithAngleAndTemplate'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("ckfTrackCandidatesP5"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(False)
)


process.ctfWithMaterialTracksP5 = cms.EDProducer("CosmicTrackSelector",
    beamspot = cms.InputTag("offlineBeamSpot"),
    chi2n_par = cms.double(10.0),
    copyExtras = cms.untracked.bool(True),
    copyTrajectories = cms.untracked.bool(False),
    keepAllTracks = cms.bool(False),
    maxNumberLostLayers = cms.uint32(999),
    max_d0 = cms.double(110.0),
    max_eta = cms.double(2.0),
    max_z0 = cms.double(300.0),
    minNumber3DLayers = cms.uint32(0),
    minNumberLayers = cms.uint32(0),
    min_nHit = cms.uint32(5),
    min_nPixelHit = cms.uint32(0),
    min_pt = cms.double(1.0),
    qualityBit = cms.string(''),
    src = cms.InputTag("ctfWithMaterialTracksCosmics")
)


process.ctfWithMaterialTracksP5Bottom = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('ctf'),
    Fitter = cms.string('FittingSmootherRKP5'),
    GeometricInnerState = cms.bool(True),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    Propagator = cms.string('RungeKuttaTrackerPropagator'),
    SimpleMagneticField = cms.string(''),
    TTRHBuilder = cms.string('WithAngleAndTemplate'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag("topBottomClusterInfoProducerBottom"),
    src = cms.InputTag("ckfTrackCandidatesP5Bottom"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(False)
)


process.ctfWithMaterialTracksP5LHCNavigation = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('ctf'),
    Fitter = cms.string('FittingSmootherRKP5'),
    GeometricInnerState = cms.bool(True),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    Propagator = cms.string('RungeKuttaTrackerPropagator'),
    SimpleMagneticField = cms.string(''),
    TTRHBuilder = cms.string('WithAngleAndTemplate'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("ckfTrackCandidatesP5LHCNavigation"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(False)
)


process.ctfWithMaterialTracksP5Top = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('ctf'),
    Fitter = cms.string('FittingSmootherRKP5'),
    GeometricInnerState = cms.bool(True),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    Propagator = cms.string('RungeKuttaTrackerPropagator'),
    SimpleMagneticField = cms.string(''),
    TTRHBuilder = cms.string('WithAngleAndTemplate'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag("topBottomClusterInfoProducerTop"),
    src = cms.InputTag("ckfTrackCandidatesP5Top"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(False)
)


process.dedxDiscrimASmi = cms.EDProducer("DeDxEstimatorProducer",
    MeVperADCPixel = cms.double(3.61e-06),
    MeVperADCStrip = cms.double(0.0009566500000000001),
    ProbabilityMode = cms.string('Accumulation'),
    Record = cms.string('SiStripDeDxMip_3D_Rcd'),
    ShapeTest = cms.bool(True),
    UseCalibration = cms.bool(False),
    UsePixel = cms.bool(False),
    UseStrip = cms.bool(True),
    calibrationPath = cms.string(''),
    estimator = cms.string('asmirnovDiscrim'),
    exponent = cms.double(-2),
    fraction = cms.double(0.4),
    mightGet = cms.optional.untracked.vstring,
    tracks = cms.InputTag("generalTracks"),
    truncate = cms.bool(True)
)


process.dedxDiscrimASmiCTF = cms.EDProducer("DeDxEstimatorProducer",
    MeVperADCPixel = cms.double(3.61e-06),
    MeVperADCStrip = cms.double(0.0009566500000000001),
    ProbabilityMode = cms.string('Accumulation'),
    Record = cms.string('SiStripDeDxMip_3D_Rcd'),
    ShapeTest = cms.bool(True),
    UseCalibration = cms.bool(False),
    UsePixel = cms.bool(False),
    UseStrip = cms.bool(True),
    calibrationPath = cms.string(''),
    estimator = cms.string('asmirnovDiscrim'),
    exponent = cms.double(-2),
    fraction = cms.double(0.4),
    mightGet = cms.optional.untracked.vstring,
    tracks = cms.InputTag("ctfWithMaterialTracksP5"),
    truncate = cms.bool(True)
)


process.dedxDiscrimASmiCTFP5LHC = cms.EDProducer("DeDxEstimatorProducer",
    MeVperADCPixel = cms.double(3.61e-06),
    MeVperADCStrip = cms.double(0.0009566500000000001),
    ProbabilityMode = cms.string('Accumulation'),
    Record = cms.string('SiStripDeDxMip_3D_Rcd'),
    ShapeTest = cms.bool(True),
    UseCalibration = cms.bool(False),
    UsePixel = cms.bool(False),
    UseStrip = cms.bool(True),
    calibrationPath = cms.string(''),
    estimator = cms.string('asmirnovDiscrim'),
    exponent = cms.double(-2),
    fraction = cms.double(0.4),
    mightGet = cms.optional.untracked.vstring,
    tracks = cms.InputTag("ctfWithMaterialTracksP5LHCNavigation"),
    truncate = cms.bool(True)
)


process.dedxDiscrimASmiCosmicTF = cms.EDProducer("DeDxEstimatorProducer",
    MeVperADCPixel = cms.double(3.61e-06),
    MeVperADCStrip = cms.double(0.0009566500000000001),
    ProbabilityMode = cms.string('Accumulation'),
    Record = cms.string('SiStripDeDxMip_3D_Rcd'),
    ShapeTest = cms.bool(True),
    UseCalibration = cms.bool(False),
    UsePixel = cms.bool(False),
    UseStrip = cms.bool(True),
    calibrationPath = cms.string(''),
    estimator = cms.string('asmirnovDiscrim'),
    exponent = cms.double(-2),
    fraction = cms.double(0.4),
    mightGet = cms.optional.untracked.vstring,
    tracks = cms.InputTag("cosmictrackfinderP5"),
    truncate = cms.bool(True)
)


process.dedxDiscrimBTag = cms.EDProducer("DeDxEstimatorProducer",
    MeVperADCPixel = cms.double(3.61e-06),
    MeVperADCStrip = cms.double(0.0009566500000000001),
    ProbabilityMode = cms.string('Accumulation'),
    Record = cms.string('SiStripDeDxMip_3D_Rcd'),
    ShapeTest = cms.bool(True),
    UseCalibration = cms.bool(False),
    UsePixel = cms.bool(False),
    UseStrip = cms.bool(True),
    calibrationPath = cms.string(''),
    estimator = cms.string('btagDiscrim'),
    exponent = cms.double(-2),
    fraction = cms.double(0.4),
    mightGet = cms.optional.untracked.vstring,
    tracks = cms.InputTag("generalTracks"),
    truncate = cms.bool(True)
)


process.dedxDiscrimProd = cms.EDProducer("DeDxEstimatorProducer",
    MeVperADCPixel = cms.double(3.61e-06),
    MeVperADCStrip = cms.double(0.0009566500000000001),
    ProbabilityMode = cms.string('Accumulation'),
    Record = cms.string('SiStripDeDxMip_3D_Rcd'),
    ShapeTest = cms.bool(True),
    UseCalibration = cms.bool(False),
    UsePixel = cms.bool(False),
    UseStrip = cms.bool(True),
    calibrationPath = cms.string(''),
    estimator = cms.string('productDiscrim'),
    exponent = cms.double(-2),
    fraction = cms.double(0.4),
    mightGet = cms.optional.untracked.vstring,
    tracks = cms.InputTag("generalTracks"),
    truncate = cms.bool(True)
)


process.dedxDiscrimSmi = cms.EDProducer("DeDxEstimatorProducer",
    MeVperADCPixel = cms.double(3.61e-06),
    MeVperADCStrip = cms.double(0.0009566500000000001),
    ProbabilityMode = cms.string('Accumulation'),
    Record = cms.string('SiStripDeDxMip_3D_Rcd'),
    ShapeTest = cms.bool(True),
    UseCalibration = cms.bool(False),
    UsePixel = cms.bool(False),
    UseStrip = cms.bool(True),
    calibrationPath = cms.string(''),
    estimator = cms.string('smirnovDiscrim'),
    exponent = cms.double(-2),
    fraction = cms.double(0.4),
    mightGet = cms.optional.untracked.vstring,
    tracks = cms.InputTag("generalTracks"),
    truncate = cms.bool(True)
)


process.dedxHarmonic2 = cms.EDProducer("DeDxEstimatorProducer",
    MeVperADCPixel = cms.double(3.61e-06),
    MeVperADCStrip = cms.double(0.0009566500000000001),
    ProbabilityMode = cms.string('Accumulation'),
    Record = cms.string('SiStripDeDxMip_3D_Rcd'),
    ShapeTest = cms.bool(True),
    UseCalibration = cms.bool(False),
    UsePixel = cms.bool(False),
    UseStrip = cms.bool(True),
    calibrationPath = cms.string(''),
    estimator = cms.string('generic'),
    exponent = cms.double(-2),
    fraction = cms.double(0.4),
    mightGet = cms.optional.untracked.vstring,
    tracks = cms.InputTag("generalTracks"),
    truncate = cms.bool(True)
)


process.dedxHarmonic2CTF = cms.EDProducer("DeDxEstimatorProducer",
    MeVperADCPixel = cms.double(3.61e-06),
    MeVperADCStrip = cms.double(0.0009566500000000001),
    ProbabilityMode = cms.string('Accumulation'),
    Record = cms.string('SiStripDeDxMip_3D_Rcd'),
    ShapeTest = cms.bool(True),
    UseCalibration = cms.bool(False),
    UsePixel = cms.bool(False),
    UseStrip = cms.bool(True),
    calibrationPath = cms.string(''),
    estimator = cms.string('generic'),
    exponent = cms.double(-2),
    fraction = cms.double(0.4),
    mightGet = cms.optional.untracked.vstring,
    tracks = cms.InputTag("ctfWithMaterialTracksP5"),
    truncate = cms.bool(True)
)


process.dedxHarmonic2CTFP5LHC = cms.EDProducer("DeDxEstimatorProducer",
    MeVperADCPixel = cms.double(3.61e-06),
    MeVperADCStrip = cms.double(0.0009566500000000001),
    ProbabilityMode = cms.string('Accumulation'),
    Record = cms.string('SiStripDeDxMip_3D_Rcd'),
    ShapeTest = cms.bool(True),
    UseCalibration = cms.bool(False),
    UsePixel = cms.bool(False),
    UseStrip = cms.bool(True),
    calibrationPath = cms.string(''),
    estimator = cms.string('generic'),
    exponent = cms.double(-2),
    fraction = cms.double(0.4),
    mightGet = cms.optional.untracked.vstring,
    tracks = cms.InputTag("ctfWithMaterialTracksP5LHCNavigation"),
    truncate = cms.bool(True)
)


process.dedxHarmonic2CosmicTF = cms.EDProducer("DeDxEstimatorProducer",
    MeVperADCPixel = cms.double(3.61e-06),
    MeVperADCStrip = cms.double(0.0009566500000000001),
    ProbabilityMode = cms.string('Accumulation'),
    Record = cms.string('SiStripDeDxMip_3D_Rcd'),
    ShapeTest = cms.bool(True),
    UseCalibration = cms.bool(False),
    UsePixel = cms.bool(False),
    UseStrip = cms.bool(True),
    calibrationPath = cms.string(''),
    estimator = cms.string('generic'),
    exponent = cms.double(-2),
    fraction = cms.double(0.4),
    mightGet = cms.optional.untracked.vstring,
    tracks = cms.InputTag("cosmictrackfinderP5"),
    truncate = cms.bool(True)
)


process.dedxHitInfo = cms.EDProducer("DeDxHitInfoProducer",
    MeVperADCPixel = cms.double(3.61e-06),
    MeVperADCStrip = cms.double(0.0009566500000000001),
    calibrationPath = cms.string('file:Gains.root'),
    lowPtTracksDeDxThreshold = cms.double(3.5),
    lowPtTracksEstimatorParameters = cms.PSet(
        exponent = cms.double(-2.0),
        fraction = cms.double(0.0),
        truncate = cms.bool(False)
    ),
    lowPtTracksPrescaleFail = cms.uint32(2000),
    lowPtTracksPrescalePass = cms.uint32(100),
    maxTrackEta = cms.double(5.0),
    minTrackHits = cms.uint32(0),
    minTrackPt = cms.double(10),
    minTrackPtPrescale = cms.double(0.5),
    shapeTest = cms.bool(True),
    tracks = cms.InputTag("generalTracks"),
    useCalibration = cms.bool(False),
    usePixel = cms.bool(True),
    usePixelForPrescales = cms.bool(False),
    useStrip = cms.bool(True)
)


process.dedxHitInfoCTF = cms.EDProducer("DeDxHitInfoProducer",
    MeVperADCPixel = cms.double(3.61e-06),
    MeVperADCStrip = cms.double(0.0009566500000000001),
    calibrationPath = cms.string('file:Gains.root'),
    lowPtTracksDeDxThreshold = cms.double(3.5),
    lowPtTracksEstimatorParameters = cms.PSet(
        exponent = cms.double(-2.0),
        fraction = cms.double(0.0),
        truncate = cms.bool(False)
    ),
    lowPtTracksPrescaleFail = cms.uint32(2000),
    lowPtTracksPrescalePass = cms.uint32(100),
    maxTrackEta = cms.double(5.0),
    minTrackHits = cms.uint32(0),
    minTrackPt = cms.double(10),
    minTrackPtPrescale = cms.double(0.5),
    shapeTest = cms.bool(True),
    tracks = cms.InputTag("ctfWithMaterialTracksP5"),
    useCalibration = cms.bool(False),
    usePixel = cms.bool(True),
    usePixelForPrescales = cms.bool(False),
    useStrip = cms.bool(True)
)


process.dedxHitInfoCTFP5LHC = cms.EDProducer("DeDxHitInfoProducer",
    MeVperADCPixel = cms.double(3.61e-06),
    MeVperADCStrip = cms.double(0.0009566500000000001),
    calibrationPath = cms.string('file:Gains.root'),
    lowPtTracksDeDxThreshold = cms.double(3.5),
    lowPtTracksEstimatorParameters = cms.PSet(
        exponent = cms.double(-2.0),
        fraction = cms.double(0.0),
        truncate = cms.bool(False)
    ),
    lowPtTracksPrescaleFail = cms.uint32(2000),
    lowPtTracksPrescalePass = cms.uint32(100),
    maxTrackEta = cms.double(5.0),
    minTrackHits = cms.uint32(0),
    minTrackPt = cms.double(10),
    minTrackPtPrescale = cms.double(0.5),
    shapeTest = cms.bool(True),
    tracks = cms.InputTag("ctfWithMaterialTracksP5LHCNavigation"),
    useCalibration = cms.bool(False),
    usePixel = cms.bool(True),
    usePixelForPrescales = cms.bool(False),
    useStrip = cms.bool(True)
)


process.dedxHitInfoCosmicTF = cms.EDProducer("DeDxHitInfoProducer",
    MeVperADCPixel = cms.double(3.61e-06),
    MeVperADCStrip = cms.double(0.0009566500000000001),
    calibrationPath = cms.string('file:Gains.root'),
    lowPtTracksDeDxThreshold = cms.double(3.5),
    lowPtTracksEstimatorParameters = cms.PSet(
        exponent = cms.double(-2.0),
        fraction = cms.double(0.0),
        truncate = cms.bool(False)
    ),
    lowPtTracksPrescaleFail = cms.uint32(2000),
    lowPtTracksPrescalePass = cms.uint32(100),
    maxTrackEta = cms.double(5.0),
    minTrackHits = cms.uint32(0),
    minTrackPt = cms.double(10),
    minTrackPtPrescale = cms.double(0.5),
    shapeTest = cms.bool(True),
    tracks = cms.InputTag("cosmictrackfinderP5"),
    useCalibration = cms.bool(False),
    usePixel = cms.bool(True),
    usePixelForPrescales = cms.bool(False),
    useStrip = cms.bool(True)
)


process.dedxMedian = cms.EDProducer("DeDxEstimatorProducer",
    MeVperADCPixel = cms.double(3.61e-06),
    MeVperADCStrip = cms.double(0.0009566500000000001),
    ProbabilityMode = cms.string('Accumulation'),
    Record = cms.string('SiStripDeDxMip_3D_Rcd'),
    ShapeTest = cms.bool(True),
    UseCalibration = cms.bool(False),
    UsePixel = cms.bool(False),
    UseStrip = cms.bool(True),
    calibrationPath = cms.string(''),
    estimator = cms.string('median'),
    exponent = cms.double(-2),
    fraction = cms.double(0.4),
    mightGet = cms.optional.untracked.vstring,
    tracks = cms.InputTag("generalTracks"),
    truncate = cms.bool(True)
)


process.dedxPixelAndStripHarmonic2T085 = cms.EDProducer("DeDxEstimatorProducer",
    MeVperADCPixel = cms.double(3.61e-06),
    MeVperADCStrip = cms.double(0.0009566500000000001),
    ProbabilityMode = cms.string('Accumulation'),
    Record = cms.string('SiStripDeDxMip_3D_Rcd'),
    ShapeTest = cms.bool(True),
    UseCalibration = cms.bool(False),
    UsePixel = cms.bool(True),
    UseStrip = cms.bool(True),
    calibrationPath = cms.string(''),
    estimator = cms.string('genericTruncated'),
    exponent = cms.double(-2),
    fraction = cms.double(-0.15),
    mightGet = cms.optional.untracked.vstring,
    tracks = cms.InputTag("generalTracks"),
    truncate = cms.bool(True)
)


process.dedxPixelHarmonic2 = cms.EDProducer("DeDxEstimatorProducer",
    MeVperADCPixel = cms.double(3.61e-06),
    MeVperADCStrip = cms.double(0.0009566500000000001),
    ProbabilityMode = cms.string('Accumulation'),
    Record = cms.string('SiStripDeDxMip_3D_Rcd'),
    ShapeTest = cms.bool(True),
    UseCalibration = cms.bool(False),
    UsePixel = cms.bool(True),
    UseStrip = cms.bool(False),
    calibrationPath = cms.string(''),
    estimator = cms.string('generic'),
    exponent = cms.double(-2),
    fraction = cms.double(0.4),
    mightGet = cms.optional.untracked.vstring,
    tracks = cms.InputTag("generalTracks"),
    truncate = cms.bool(True)
)


process.dedxTruncated40 = cms.EDProducer("DeDxEstimatorProducer",
    MeVperADCPixel = cms.double(3.61e-06),
    MeVperADCStrip = cms.double(0.0009566500000000001),
    ProbabilityMode = cms.string('Accumulation'),
    Record = cms.string('SiStripDeDxMip_3D_Rcd'),
    ShapeTest = cms.bool(True),
    UseCalibration = cms.bool(False),
    UsePixel = cms.bool(False),
    UseStrip = cms.bool(True),
    calibrationPath = cms.string(''),
    estimator = cms.string('truncated'),
    exponent = cms.double(-2),
    fraction = cms.double(0.4),
    mightGet = cms.optional.untracked.vstring,
    tracks = cms.InputTag("generalTracks"),
    truncate = cms.bool(True)
)


process.dedxTruncated40CTF = cms.EDProducer("DeDxEstimatorProducer",
    MeVperADCPixel = cms.double(3.61e-06),
    MeVperADCStrip = cms.double(0.0009566500000000001),
    ProbabilityMode = cms.string('Accumulation'),
    Record = cms.string('SiStripDeDxMip_3D_Rcd'),
    ShapeTest = cms.bool(True),
    UseCalibration = cms.bool(False),
    UsePixel = cms.bool(False),
    UseStrip = cms.bool(True),
    calibrationPath = cms.string(''),
    estimator = cms.string('truncated'),
    exponent = cms.double(-2),
    fraction = cms.double(0.4),
    mightGet = cms.optional.untracked.vstring,
    tracks = cms.InputTag("ctfWithMaterialTracksP5"),
    truncate = cms.bool(True)
)


process.dedxTruncated40CTFP5LHC = cms.EDProducer("DeDxEstimatorProducer",
    MeVperADCPixel = cms.double(3.61e-06),
    MeVperADCStrip = cms.double(0.0009566500000000001),
    ProbabilityMode = cms.string('Accumulation'),
    Record = cms.string('SiStripDeDxMip_3D_Rcd'),
    ShapeTest = cms.bool(True),
    UseCalibration = cms.bool(False),
    UsePixel = cms.bool(False),
    UseStrip = cms.bool(True),
    calibrationPath = cms.string(''),
    estimator = cms.string('truncated'),
    exponent = cms.double(-2),
    fraction = cms.double(0.4),
    mightGet = cms.optional.untracked.vstring,
    tracks = cms.InputTag("ctfWithMaterialTracksP5LHCNavigation"),
    truncate = cms.bool(True)
)


process.dedxTruncated40CosmicTF = cms.EDProducer("DeDxEstimatorProducer",
    MeVperADCPixel = cms.double(3.61e-06),
    MeVperADCStrip = cms.double(0.0009566500000000001),
    ProbabilityMode = cms.string('Accumulation'),
    Record = cms.string('SiStripDeDxMip_3D_Rcd'),
    ShapeTest = cms.bool(True),
    UseCalibration = cms.bool(False),
    UsePixel = cms.bool(False),
    UseStrip = cms.bool(True),
    calibrationPath = cms.string(''),
    estimator = cms.string('truncated'),
    exponent = cms.double(-2),
    fraction = cms.double(0.4),
    mightGet = cms.optional.untracked.vstring,
    tracks = cms.InputTag("cosmictrackfinderP5"),
    truncate = cms.bool(True)
)


process.dedxUnbinned = cms.EDProducer("DeDxEstimatorProducer",
    MeVperADCPixel = cms.double(3.61e-06),
    MeVperADCStrip = cms.double(0.0009566500000000001),
    ProbabilityMode = cms.string('Accumulation'),
    Record = cms.string('SiStripDeDxMip_3D_Rcd'),
    ShapeTest = cms.bool(True),
    UseCalibration = cms.bool(False),
    UsePixel = cms.bool(False),
    UseStrip = cms.bool(True),
    calibrationPath = cms.string(''),
    estimator = cms.string('unbinnedFit'),
    exponent = cms.double(-2),
    fraction = cms.double(0.4),
    mightGet = cms.optional.untracked.vstring,
    tracks = cms.InputTag("generalTracks"),
    truncate = cms.bool(True)
)


process.detachedQuadStep = cms.EDProducer("TrackTfClassifier",
    beamspot = cms.InputTag("offlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mightGet = cms.optional.untracked.vstring,
    mva = cms.PSet(
        batchSize = cms.int32(16),
        tfDnnLabel = cms.string('trackSelectionTf')
    ),
    qualityCuts = cms.vdouble(-0.66, -0.15, 0.46),
    src = cms.InputTag("detachedQuadStepTracks"),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.detachedQuadStepClusters = cms.EDProducer("TrackClusterRemover",
    TrackQuality = cms.string('highPurity'),
    maxChi2 = cms.double(9.0),
    mightGet = cms.optional.untracked.vstring,
    minNumberOfLayersWithMeasBeforeFiltering = cms.int32(0),
    oldClusterRemovalInfo = cms.InputTag("lowPtTripletStepClusters"),
    overrideTrkQuals = cms.InputTag(""),
    pixelClusters = cms.InputTag("siPixelClusters"),
    stripClusters = cms.InputTag("siStripClusters"),
    trackClassifier = cms.InputTag("lowPtTripletStep","QualityMasks"),
    trajectories = cms.InputTag("lowPtTripletStepTracks")
)


process.detachedQuadStepHitDoublets = cms.EDProducer("HitPairEDProducer",
    clusterCheck = cms.InputTag("trackerClusterCheck"),
    layerPairs = cms.vuint32(0, 1, 2),
    maxElement = cms.uint32(50000000),
    maxElementTotal = cms.uint32(50000000),
    mightGet = cms.optional.untracked.vstring,
    produceIntermediateHitDoublets = cms.bool(True),
    produceSeedingHitSets = cms.bool(False),
    seedingLayers = cms.InputTag("detachedQuadStepSeedLayers"),
    trackingRegions = cms.InputTag("detachedQuadStepTrackingRegions"),
    trackingRegionsSeedingLayers = cms.InputTag("")
)


process.detachedQuadStepHitQuadruplets = cms.EDProducer("CAHitQuadrupletEDProducer",
    CAHardPtCut = cms.double(0),
    CAOnlyOneLastHitPerLayerFilter = cms.optional.bool,
    CAPhiCut = cms.double(0),
    CAPhiCut_byTriplets = cms.VPSet(cms.PSet(
        cut = cms.double(-1),
        seedingLayers = cms.string('')
    )),
    CAThetaCut = cms.double(0.0011),
    CAThetaCut_byTriplets = cms.VPSet(cms.PSet(
        cut = cms.double(-1),
        seedingLayers = cms.string('')
    )),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('none')
    ),
    doublets = cms.InputTag("detachedQuadStepHitDoublets"),
    extraHitRPhitolerance = cms.double(0),
    fitFastCircle = cms.bool(True),
    fitFastCircleChi2Cut = cms.bool(True),
    maxChi2 = cms.PSet(
        enabled = cms.bool(True),
        pt1 = cms.double(0.8),
        pt2 = cms.double(2),
        value1 = cms.double(500),
        value2 = cms.double(100)
    ),
    mightGet = cms.optional.untracked.vstring,
    useBendingCorrection = cms.bool(True)
)


process.detachedQuadStepMasks = cms.EDProducer("FastTrackerRecHitMaskProducer",
    TrackQuality = cms.string('highPurity'),
    minNumberOfLayersWithMeasBeforeFiltering = cms.int32(0),
    minNumberOfLayersWithMeasBeforeFiltering_ = cms.int32(0),
    oldHitRemovalInfo = cms.InputTag("lowPtTripletStepMasks"),
    recHits = cms.InputTag("fastTrackerRecHits"),
    trackClassifier = cms.InputTag("lowPtTripletStep","QualityMasks"),
    trackQuality = cms.string('hightPurity'),
    trajectories = cms.InputTag("lowPtTripletStepTracks")
)


process.detachedQuadStepSeedLayers = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(
        HitProducer = cms.string('siPixelRecHits'),
        TTRHBuilder = cms.string('WithTrackAngle'),
        skipClusters = cms.InputTag("detachedQuadStepClusters")
    ),
    FPix = cms.PSet(
        HitProducer = cms.string('siPixelRecHits'),
        TTRHBuilder = cms.string('WithTrackAngle'),
        skipClusters = cms.InputTag("detachedQuadStepClusters")
    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(

    ),
    TIB = cms.PSet(

    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(

    ),
    layerList = cms.vstring(
        'BPix1+BPix2+BPix3+BPix4',
        'BPix1+BPix2+BPix3+FPix1_pos',
        'BPix1+BPix2+BPix3+FPix1_neg',
        'BPix1+BPix2+FPix1_pos+FPix2_pos',
        'BPix1+BPix2+FPix1_neg+FPix2_neg',
        'BPix1+FPix1_pos+FPix2_pos+FPix3_pos',
        'BPix1+FPix1_neg+FPix2_neg+FPix3_neg'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.detachedQuadStepSeeds = cms.EDProducer("SeedCreatorFromRegionConsecutiveHitsTripletOnlyEDProducer",
    MinOneOverPtError = cms.double(1),
    OriginTransverseErrorMultiplier = cms.double(1),
    SeedComparitorPSet = cms.PSet(
        ClusterShapeCacheSrc = cms.InputTag("siPixelClusterShapeCache"),
        ClusterShapeHitFilterName = cms.string('ClusterShapeHitFilter'),
        ComponentName = cms.string('PixelClusterShapeSeedComparitor'),
        FilterAtHelixStage = cms.bool(False),
        FilterPixelHits = cms.bool(True),
        FilterStripHits = cms.bool(False)
    ),
    SeedMomentumForBOFF = cms.double(5),
    TTRHBuilder = cms.string('WithTrackAngle'),
    forceKinematicWithRegionDirection = cms.bool(False),
    magneticField = cms.string('ParabolicMf'),
    mightGet = cms.optional.untracked.vstring,
    propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    seedingHitSets = cms.InputTag("detachedQuadStepHitQuadruplets")
)


process.detachedQuadStepSelector = cms.EDProducer("MultiTrackSelector",
    beamspot = cms.InputTag("offlineBeamSpot"),
    src = cms.InputTag("detachedQuadStepTracks"),
    trackSelectors = cms.VPSet(
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(True),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(1.0),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(0.9, 3.0),
            d0_par2 = cms.vdouble(1.0, 3.0),
            dz_par1 = cms.vdouble(0.9, 3.0),
            dz_par2 = cms.vdouble(1.0, 3.0),
            keepAllTracks = cms.bool(False),
            maxNumberLostLayers = cms.uint32(999),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(0),
            minNumberLayers = cms.uint32(3),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(0),
            nSigmaZ = cms.double(4.0),
            name = cms.string('detachedQuadStepVtxLoose'),
            preFilterName = cms.string(''),
            qualityBit = cms.string('loose'),
            res_par = cms.vdouble(0.003, 0.001),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        ),
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(True),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(0.6),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(1.3, 4.0),
            d0_par2 = cms.vdouble(1.3, 4.0),
            dz_par1 = cms.vdouble(1.3, 4.0),
            dz_par2 = cms.vdouble(1.3, 4.0),
            keepAllTracks = cms.bool(False),
            maxNumberLostLayers = cms.uint32(999),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(0),
            minNumberLayers = cms.uint32(3),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(0),
            nSigmaZ = cms.double(4.0),
            name = cms.string('detachedQuadStepTrkLoose'),
            preFilterName = cms.string(''),
            qualityBit = cms.string('loose'),
            res_par = cms.vdouble(0.003, 0.001),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        ),
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(True),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(0.9),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(0.9, 3.0),
            d0_par2 = cms.vdouble(0.9, 3.0),
            dz_par1 = cms.vdouble(0.9, 3.0),
            dz_par2 = cms.vdouble(0.9, 3.0),
            keepAllTracks = cms.bool(True),
            maxNumberLostLayers = cms.uint32(1),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(3),
            minNumberLayers = cms.uint32(3),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(0),
            nSigmaZ = cms.double(4.0),
            name = cms.string('detachedQuadStepVtxTight'),
            preFilterName = cms.string('detachedQuadStepVtxLoose'),
            qualityBit = cms.string('tight'),
            res_par = cms.vdouble(0.003, 0.001),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        ),
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(True),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(0.5),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(1.1, 4.0),
            d0_par2 = cms.vdouble(1.1, 4.0),
            dz_par1 = cms.vdouble(1.1, 4.0),
            dz_par2 = cms.vdouble(1.1, 4.0),
            keepAllTracks = cms.bool(True),
            maxNumberLostLayers = cms.uint32(1),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(3),
            minNumberLayers = cms.uint32(4),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(0),
            nSigmaZ = cms.double(4.0),
            name = cms.string('detachedQuadStepTrkTight'),
            preFilterName = cms.string('detachedQuadStepTrkLoose'),
            qualityBit = cms.string('tight'),
            res_par = cms.vdouble(0.003, 0.001),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        ),
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(True),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(0.9),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(0.8, 3.0),
            d0_par2 = cms.vdouble(0.8, 3.0),
            dz_par1 = cms.vdouble(0.8, 3.0),
            dz_par2 = cms.vdouble(0.8, 3.0),
            keepAllTracks = cms.bool(True),
            maxNumberLostLayers = cms.uint32(1),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(4.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(3),
            minNumberLayers = cms.uint32(3),
            min_eta = cms.double(-4.0),
            min_nhits = cms.uint32(0),
            nSigmaZ = cms.double(4.0),
            name = cms.string('detachedQuadStepVtx'),
            preFilterName = cms.string('detachedQuadStepVtxTight'),
            qualityBit = cms.string('highPurity'),
            res_par = cms.vdouble(0.003, 0.001),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        ),
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(True),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(0.5),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(0.9, 4.0),
            d0_par2 = cms.vdouble(0.9, 4.0),
            dz_par1 = cms.vdouble(0.9, 4.0),
            dz_par2 = cms.vdouble(0.9, 4.0),
            keepAllTracks = cms.bool(True),
            maxNumberLostLayers = cms.uint32(1),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(4.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(3),
            minNumberLayers = cms.uint32(4),
            min_eta = cms.double(-4.0),
            min_nhits = cms.uint32(0),
            nSigmaZ = cms.double(4.0),
            name = cms.string('detachedQuadStepTrk'),
            preFilterName = cms.string('detachedQuadStepTrkTight'),
            qualityBit = cms.string('highPurity'),
            res_par = cms.vdouble(0.003, 0.001),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        )
    ),
    useVertices = cms.bool(True),
    useVtxError = cms.bool(False),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.detachedQuadStepTrackCandidates = cms.EDProducer("MkFitOutputConverter",
    batchSize = cms.int32(16),
    candMVASel = cms.bool(False),
    candWP = cms.double(0),
    doErrorRescale = cms.bool(True),
    mightGet = cms.optional.untracked.vstring,
    mkFitEventOfHits = cms.InputTag("mkFitEventOfHits"),
    mkFitPixelHits = cms.InputTag("mkFitSiPixelHits"),
    mkFitSeeds = cms.InputTag("detachedQuadStepTrackCandidatesMkFitSeeds"),
    mkFitStripHits = cms.InputTag("mkFitSiStripHits"),
    propagatorAlong = cms.ESInputTag("","PropagatorWithMaterial"),
    propagatorOpposite = cms.ESInputTag("","PropagatorWithMaterialOpposite"),
    qualityMaxInvPt = cms.double(100),
    qualityMaxPosErr = cms.double(100),
    qualityMaxR = cms.double(120),
    qualityMaxZ = cms.double(280),
    qualityMinTheta = cms.double(0.01),
    qualitySignPt = cms.bool(True),
    seeds = cms.InputTag("detachedQuadStepSeeds"),
    tfDnnLabel = cms.string('trackSelectionTf'),
    tracks = cms.InputTag("detachedQuadStepTrackCandidatesMkFit"),
    ttrhBuilder = cms.ESInputTag("","WithTrackAngle")
)


process.detachedQuadStepTrackCandidatesMkFit = cms.EDProducer("MkFitProducer",
    backwardFitInCMSSW = cms.bool(False),
    buildingRoutine = cms.string('cloneEngine'),
    clustersToSkip = cms.InputTag("detachedQuadStepClusters"),
    config = cms.ESInputTag("","detachedQuadStepTrackCandidatesMkFitConfig"),
    eventOfHits = cms.InputTag("mkFitEventOfHits"),
    limitConcurrency = cms.untracked.bool(False),
    mightGet = cms.optional.untracked.vstring,
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutLoose')
    ),
    mkFitSilent = cms.untracked.bool(True),
    pixelHits = cms.InputTag("mkFitSiPixelHits"),
    removeDuplicates = cms.bool(True),
    seedCleaning = cms.bool(True),
    seeds = cms.InputTag("detachedQuadStepTrackCandidatesMkFitSeeds"),
    stripHits = cms.InputTag("mkFitSiStripHits")
)


process.detachedQuadStepTrackCandidatesMkFitSeeds = cms.EDProducer("MkFitSeedConverter",
    maxNSeeds = cms.uint32(500000),
    mightGet = cms.optional.untracked.vstring,
    seeds = cms.InputTag("detachedQuadStepSeeds"),
    ttrhBuilder = cms.ESInputTag("","WithTrackAngle")
)


process.detachedQuadStepTrackingRegions = cms.EDProducer("GlobalTrackingRegionFromBeamSpotEDProducer",
    RegionPSet = cms.PSet(
        beamSpot = cms.InputTag("offlineBeamSpot"),
        nSigmaZ = cms.double(0),
        originHalfLength = cms.double(15.0),
        originRadius = cms.double(1.5),
        precise = cms.bool(True),
        ptMin = cms.double(0.3),
        useMultipleScattering = cms.bool(False)
    ),
    mightGet = cms.optional.untracked.vstring
)


process.detachedQuadStepTracks = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('detachedQuadStep'),
    Fitter = cms.string('FlexibleKFFittingSmoother'),
    GeometricInnerState = cms.bool(False),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    Propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    TTRHBuilder = cms.string('WithAngleAndTemplateWithoutProbQ'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("detachedQuadStepTrackCandidates"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(True)
)


process.detachedTripletStep = cms.EDProducer("TrackTfClassifier",
    beamspot = cms.InputTag("offlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mightGet = cms.optional.untracked.vstring,
    mva = cms.PSet(
        batchSize = cms.int32(16),
        tfDnnLabel = cms.string('trackSelectionTf')
    ),
    qualityCuts = cms.vdouble(-0.42, 0.16, 0.78),
    src = cms.InputTag("detachedTripletStepTracks"),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.detachedTripletStepClassifier1 = cms.EDProducer("TrackMVAClassifierDetached",
    beamspot = cms.InputTag("offlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mightGet = cms.optional.untracked.vstring,
    mva = cms.PSet(
        GBRForestFileName = cms.string(''),
        GBRForestLabel = cms.string('MVASelectorIter3_13TeV')
    ),
    qualityCuts = cms.vdouble(-0.5, 0.0, 0.5),
    src = cms.InputTag("detachedTripletStepTracks"),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.detachedTripletStepClassifier2 = cms.EDProducer("TrackMVAClassifierPrompt",
    beamspot = cms.InputTag("offlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mightGet = cms.optional.untracked.vstring,
    mva = cms.PSet(
        GBRForestFileName = cms.string(''),
        GBRForestLabel = cms.string('MVASelectorIter0_13TeV')
    ),
    qualityCuts = cms.vdouble(-0.2, 0.0, 0.4),
    src = cms.InputTag("detachedTripletStepTracks"),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.detachedTripletStepClusters = cms.EDProducer("TrackClusterRemover",
    TrackQuality = cms.string('highPurity'),
    maxChi2 = cms.double(9.0),
    mightGet = cms.optional.untracked.vstring,
    minNumberOfLayersWithMeasBeforeFiltering = cms.int32(0),
    oldClusterRemovalInfo = cms.InputTag("detachedQuadStepClusters"),
    overrideTrkQuals = cms.InputTag(""),
    pixelClusters = cms.InputTag("siPixelClusters"),
    stripClusters = cms.InputTag("siStripClusters"),
    trackClassifier = cms.InputTag("detachedQuadStep","QualityMasks"),
    trajectories = cms.InputTag("detachedQuadStepTracks")
)


process.detachedTripletStepHitDoublets = cms.EDProducer("HitPairEDProducer",
    clusterCheck = cms.InputTag("trackerClusterCheck"),
    layerPairs = cms.vuint32(0, 1),
    maxElement = cms.uint32(50000000),
    maxElementTotal = cms.uint32(50000000),
    mightGet = cms.optional.untracked.vstring,
    produceIntermediateHitDoublets = cms.bool(True),
    produceSeedingHitSets = cms.bool(False),
    seedingLayers = cms.InputTag("detachedTripletStepSeedLayers"),
    trackingRegions = cms.InputTag("detachedTripletStepTrackingRegions"),
    trackingRegionsSeedingLayers = cms.InputTag("")
)


process.detachedTripletStepHitTriplets = cms.EDProducer("CAHitTripletEDProducer",
    CAHardPtCut = cms.double(0.2),
    CAPhiCut = cms.double(0),
    CAPhiCut_byTriplets = cms.VPSet(cms.PSet(
        cut = cms.double(-1),
        seedingLayers = cms.string('')
    )),
    CAThetaCut = cms.double(0.001),
    CAThetaCut_byTriplets = cms.VPSet(cms.PSet(
        cut = cms.double(-1),
        seedingLayers = cms.string('')
    )),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('none')
    ),
    doublets = cms.InputTag("detachedTripletStepHitDoublets"),
    extraHitRPhitolerance = cms.double(0),
    maxChi2 = cms.PSet(
        enabled = cms.bool(True),
        pt1 = cms.double(0.8),
        pt2 = cms.double(2),
        value1 = cms.double(300),
        value2 = cms.double(10)
    ),
    mightGet = cms.optional.untracked.vstring,
    useBendingCorrection = cms.bool(True)
)


process.detachedTripletStepMasks = cms.EDProducer("FastTrackerRecHitMaskProducer",
    TrackQuality = cms.string('highPurity'),
    minNumberOfLayersWithMeasBeforeFiltering = cms.int32(0),
    minNumberOfLayersWithMeasBeforeFiltering_ = cms.int32(0),
    oldHitRemovalInfo = cms.InputTag("detachedQuadStepMasks"),
    recHits = cms.InputTag("fastTrackerRecHits"),
    trackClassifier = cms.InputTag("detachedQuadStep","QualityMasks"),
    trackQuality = cms.string('hightPurity'),
    trajectories = cms.InputTag("detachedQuadStepTracks")
)


process.detachedTripletStepSeedClusterMask = cms.EDProducer("SeedClusterRemover",
    Common = cms.PSet(
        maxChi2 = cms.double(9.0)
    ),
    oldClusterRemovalInfo = cms.InputTag("initialStepSeedClusterMask"),
    pixelClusters = cms.InputTag("siPixelClusters"),
    stripClusters = cms.InputTag("siStripClusters"),
    trajectories = cms.InputTag("lowPtTripletStepSeeds")
)


process.detachedTripletStepSeedLayers = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(
        HitProducer = cms.string('siPixelRecHits'),
        TTRHBuilder = cms.string('WithTrackAngle'),
        skipClusters = cms.InputTag("detachedTripletStepClusters")
    ),
    FPix = cms.PSet(
        HitProducer = cms.string('siPixelRecHits'),
        TTRHBuilder = cms.string('WithTrackAngle'),
        skipClusters = cms.InputTag("detachedTripletStepClusters")
    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(

    ),
    TIB = cms.PSet(

    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(

    ),
    layerList = cms.vstring(
        'BPix1+BPix2+BPix3',
        'BPix2+BPix3+BPix4',
        'BPix2+BPix3+FPix1_pos',
        'BPix2+BPix3+FPix1_neg',
        'BPix2+FPix1_pos+FPix2_pos',
        'BPix2+FPix1_neg+FPix2_neg',
        'FPix1_pos+FPix2_pos+FPix3_pos',
        'FPix1_neg+FPix2_neg+FPix3_neg'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.detachedTripletStepSeeds = cms.EDProducer("SeedCreatorFromRegionConsecutiveHitsTripletOnlyEDProducer",
    MinOneOverPtError = cms.double(1),
    OriginTransverseErrorMultiplier = cms.double(1),
    SeedComparitorPSet = cms.PSet(
        ClusterShapeCacheSrc = cms.InputTag("siPixelClusterShapeCache"),
        ClusterShapeHitFilterName = cms.string('ClusterShapeHitFilter'),
        ComponentName = cms.string('PixelClusterShapeSeedComparitor'),
        FilterAtHelixStage = cms.bool(False),
        FilterPixelHits = cms.bool(True),
        FilterStripHits = cms.bool(False)
    ),
    SeedMomentumForBOFF = cms.double(5),
    TTRHBuilder = cms.string('WithTrackAngle'),
    forceKinematicWithRegionDirection = cms.bool(False),
    magneticField = cms.string('ParabolicMf'),
    mightGet = cms.optional.untracked.vstring,
    propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    seedingHitSets = cms.InputTag("detachedTripletStepHitTriplets")
)


process.detachedTripletStepSelector = cms.EDProducer("MultiTrackSelector",
    GBRForestLabel = cms.string('MVASelectorIter3'),
    beamspot = cms.InputTag("offlineBeamSpot"),
    src = cms.InputTag("detachedTripletStepTracks"),
    trackSelectors = cms.VPSet(
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(True),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(1.6),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(1.2, 3.0),
            d0_par2 = cms.vdouble(1.3, 3.0),
            dz_par1 = cms.vdouble(1.2, 3.0),
            dz_par2 = cms.vdouble(1.3, 3.0),
            keepAllTracks = cms.bool(False),
            maxNumberLostLayers = cms.uint32(999),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(0),
            minNumberLayers = cms.uint32(3),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(0),
            nSigmaZ = cms.double(4.0),
            name = cms.string('detachedTripletStepVtxLoose'),
            preFilterName = cms.string(''),
            qualityBit = cms.string('loose'),
            res_par = cms.vdouble(0.003, 0.001),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        ),
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(True),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(0.7),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(1.6, 4.0),
            d0_par2 = cms.vdouble(1.6, 4.0),
            dz_par1 = cms.vdouble(1.6, 4.0),
            dz_par2 = cms.vdouble(1.6, 4.0),
            keepAllTracks = cms.bool(False),
            maxNumberLostLayers = cms.uint32(999),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(0),
            minNumberLayers = cms.uint32(3),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(0),
            nSigmaZ = cms.double(4.0),
            name = cms.string('detachedTripletStepTrkLoose'),
            preFilterName = cms.string(''),
            qualityBit = cms.string('loose'),
            res_par = cms.vdouble(0.003, 0.001),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        ),
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(True),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(0.7),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(0.95, 3.0),
            d0_par2 = cms.vdouble(1.0, 3.0),
            dz_par1 = cms.vdouble(0.9, 3.0),
            dz_par2 = cms.vdouble(1.0, 3.0),
            keepAllTracks = cms.bool(True),
            maxNumberLostLayers = cms.uint32(1),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(3),
            minNumberLayers = cms.uint32(3),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(0),
            nSigmaZ = cms.double(4.0),
            name = cms.string('detachedTripletStepVtxTight'),
            preFilterName = cms.string('detachedTripletStepVtxLoose'),
            qualityBit = cms.string('tight'),
            res_par = cms.vdouble(0.003, 0.001),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        ),
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(True),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(0.5),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(1.1, 4.0),
            d0_par2 = cms.vdouble(1.1, 4.0),
            dz_par1 = cms.vdouble(1.1, 4.0),
            dz_par2 = cms.vdouble(1.1, 4.0),
            keepAllTracks = cms.bool(True),
            maxNumberLostLayers = cms.uint32(1),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(3),
            minNumberLayers = cms.uint32(5),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(0),
            nSigmaZ = cms.double(4.0),
            name = cms.string('detachedTripletStepTrkTight'),
            preFilterName = cms.string('detachedTripletStepTrkLoose'),
            qualityBit = cms.string('tight'),
            res_par = cms.vdouble(0.003, 0.001),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        ),
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(True),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(0.7),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(0.85, 3.0),
            d0_par2 = cms.vdouble(0.9, 3.0),
            dz_par1 = cms.vdouble(0.8, 3.0),
            dz_par2 = cms.vdouble(0.9, 3.0),
            keepAllTracks = cms.bool(True),
            maxNumberLostLayers = cms.uint32(1),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(3),
            minNumberLayers = cms.uint32(3),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(0),
            nSigmaZ = cms.double(4.0),
            name = cms.string('detachedTripletStepVtx'),
            preFilterName = cms.string('detachedTripletStepVtxTight'),
            qualityBit = cms.string('highPurity'),
            res_par = cms.vdouble(0.003, 0.001),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        ),
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(True),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(0.4),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(1.0, 4.0),
            d0_par2 = cms.vdouble(1.0, 4.0),
            dz_par1 = cms.vdouble(1.0, 4.0),
            dz_par2 = cms.vdouble(1.0, 4.0),
            keepAllTracks = cms.bool(True),
            maxNumberLostLayers = cms.uint32(1),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(4),
            minNumberLayers = cms.uint32(5),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(0),
            nSigmaZ = cms.double(4.0),
            name = cms.string('detachedTripletStepTrk'),
            preFilterName = cms.string('detachedTripletStepTrkTight'),
            qualityBit = cms.string('highPurity'),
            res_par = cms.vdouble(0.003, 0.001),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        )
    ),
    useAnyMVA = cms.bool(False),
    useVertices = cms.bool(True),
    useVtxError = cms.bool(False),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.detachedTripletStepTrackCandidates = cms.EDProducer("MkFitOutputConverter",
    batchSize = cms.int32(16),
    candMVASel = cms.bool(False),
    candWP = cms.double(0),
    doErrorRescale = cms.bool(True),
    mightGet = cms.optional.untracked.vstring,
    mkFitEventOfHits = cms.InputTag("mkFitEventOfHits"),
    mkFitPixelHits = cms.InputTag("mkFitSiPixelHits"),
    mkFitSeeds = cms.InputTag("detachedTripletStepTrackCandidatesMkFitSeeds"),
    mkFitStripHits = cms.InputTag("mkFitSiStripHits"),
    propagatorAlong = cms.ESInputTag("","PropagatorWithMaterial"),
    propagatorOpposite = cms.ESInputTag("","PropagatorWithMaterialOpposite"),
    qualityMaxInvPt = cms.double(100),
    qualityMaxPosErr = cms.double(100),
    qualityMaxR = cms.double(120),
    qualityMaxZ = cms.double(280),
    qualityMinTheta = cms.double(0.01),
    qualitySignPt = cms.bool(True),
    seeds = cms.InputTag("detachedTripletStepSeeds"),
    tfDnnLabel = cms.string('trackSelectionTf'),
    tracks = cms.InputTag("detachedTripletStepTrackCandidatesMkFit"),
    ttrhBuilder = cms.ESInputTag("","WithTrackAngle")
)


process.detachedTripletStepTrackCandidatesMkFit = cms.EDProducer("MkFitProducer",
    backwardFitInCMSSW = cms.bool(False),
    buildingRoutine = cms.string('cloneEngine'),
    clustersToSkip = cms.InputTag("detachedTripletStepClusters"),
    config = cms.ESInputTag("","detachedTripletStepTrackCandidatesMkFitConfig"),
    eventOfHits = cms.InputTag("mkFitEventOfHits"),
    limitConcurrency = cms.untracked.bool(False),
    mightGet = cms.optional.untracked.vstring,
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutLoose')
    ),
    mkFitSilent = cms.untracked.bool(True),
    pixelHits = cms.InputTag("mkFitSiPixelHits"),
    removeDuplicates = cms.bool(True),
    seedCleaning = cms.bool(True),
    seeds = cms.InputTag("detachedTripletStepTrackCandidatesMkFitSeeds"),
    stripHits = cms.InputTag("mkFitSiStripHits")
)


process.detachedTripletStepTrackCandidatesMkFitSeeds = cms.EDProducer("MkFitSeedConverter",
    maxNSeeds = cms.uint32(500000),
    mightGet = cms.optional.untracked.vstring,
    seeds = cms.InputTag("detachedTripletStepSeeds"),
    ttrhBuilder = cms.ESInputTag("","WithTrackAngle")
)


process.detachedTripletStepTrackingRegions = cms.EDProducer("GlobalTrackingRegionFromBeamSpotEDProducer",
    RegionPSet = cms.PSet(
        beamSpot = cms.InputTag("offlineBeamSpot"),
        nSigmaZ = cms.double(0),
        originHalfLength = cms.double(15.0),
        originRadius = cms.double(1.5),
        precise = cms.bool(True),
        ptMin = cms.double(0.25),
        useMultipleScattering = cms.bool(False)
    ),
    mightGet = cms.optional.untracked.vstring
)


process.detachedTripletStepTracks = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('detachedTripletStep'),
    Fitter = cms.string('FlexibleKFFittingSmoother'),
    GeometricInnerState = cms.bool(False),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    Propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    TTRHBuilder = cms.string('WithAngleAndTemplateWithoutProbQ'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("detachedTripletStepTrackCandidates"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(True)
)


process.discriminationByIsolationMVArun2v1 = cms.EDProducer("RecoTauDiscriminantCutMultiplexer",
    PFTauProducer = cms.InputTag("pfTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            Producer = cms.InputTag("pfRecoTauDiscriminationByLeadingTrackFinding"),
            cut = cms.double(0.5)
        )
    ),
    inputFileName = cms.FileInPath('RecoTauTag/RecoTau/data/emptyMVAinputFile'),
    loadMVAfromDB = cms.bool(True),
    mapping = cms.VPSet(cms.PSet(
        category = cms.uint32(0),
        cut = cms.string('newDMwLT'),
        variable = cms.string('pt')
    )),
    mightGet = cms.optional.untracked.vstring,
    mvaOutput_normalization = cms.string(''),
    rawValues = cms.vstring('discriminator'),
    toMultiplex = cms.InputTag("discriminationByIsolationMVArun2v1raw"),
    verbosity = cms.int32(0),
    workingPoints = cms.vstring(
        'Eff80',
        'Eff70',
        'Eff60',
        'Eff50',
        'Eff40'
    )
)


process.discriminationByIsolationMVArun2v1raw = cms.EDProducer("PFRecoTauDiscriminationByMVAIsolationRun2",
    PFTauProducer = cms.InputTag("pfTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            Producer = cms.InputTag("pfRecoTauDiscriminationByLeadingTrackFinding"),
            cut = cms.double(0.5)
        )
    ),
    inputFileName = cms.optional.FileInPath,
    inputIDNameSuffix = cms.string(''),
    loadMVAfromDB = cms.bool(True),
    mightGet = cms.optional.untracked.vstring,
    mvaName = cms.string('tauIdMVAnewDMwLT'),
    mvaOpt = cms.string('newDMwLT'),
    srcBasicTauDiscriminators = cms.InputTag("hpsPFTauBasicDiscriminators"),
    srcTauTransverseImpactParameters = cms.InputTag(""),
    verbosity = cms.int32(0)
)


process.displacedGeneralStep = cms.EDProducer("TrackMVAClassifierDetached",
    beamspot = cms.InputTag("offlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mightGet = cms.optional.untracked.vstring,
    mva = cms.PSet(
        GBRForestFileName = cms.string(''),
        GBRForestLabel = cms.string('MVASelectorTobTecStep_Phase1')
    ),
    qualityCuts = cms.vdouble(-0.6, -0.45, -0.3),
    src = cms.InputTag("displacedGeneralStepTracks"),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.displacedGeneralStepClassifier1 = cms.EDProducer("TrackMVAClassifierDetached",
    beamspot = cms.InputTag("offlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mightGet = cms.optional.untracked.vstring,
    mva = cms.PSet(
        GBRForestFileName = cms.string(''),
        GBRForestLabel = cms.string('MVASelectorIter6_13TeV')
    ),
    qualityCuts = cms.vdouble(-0.6, -0.45, -0.3),
    src = cms.InputTag("displacedGeneralStepTracks"),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.displacedGeneralStepClassifier2 = cms.EDProducer("TrackMVAClassifierPrompt",
    beamspot = cms.InputTag("offlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mightGet = cms.optional.untracked.vstring,
    mva = cms.PSet(
        GBRForestFileName = cms.string(''),
        GBRForestLabel = cms.string('MVASelectorIter0_13TeV')
    ),
    qualityCuts = cms.vdouble(0.0, 0.0, 0.0),
    src = cms.InputTag("displacedGeneralStepTracks"),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.displacedGeneralStepClusters = cms.EDProducer("TrackClusterRemover",
    TrackQuality = cms.string('highPurity'),
    maxChi2 = cms.double(9.0),
    mightGet = cms.optional.untracked.vstring,
    minNumberOfLayersWithMeasBeforeFiltering = cms.int32(0),
    oldClusterRemovalInfo = cms.InputTag(""),
    overrideTrkQuals = cms.InputTag(""),
    pixelClusters = cms.InputTag("siPixelClusters"),
    stripClusters = cms.InputTag("siStripClusters"),
    trackClassifier = cms.InputTag("","QualityMasks"),
    trajectories = cms.InputTag("")
)


process.displacedGeneralStepHitDoublets = cms.EDProducer("HitPairEDProducer",
    clusterCheck = cms.InputTag("trackerClusterCheck"),
    layerPairs = cms.vuint32(0),
    maxElement = cms.uint32(50000000),
    maxElementTotal = cms.uint32(50000000),
    mightGet = cms.optional.untracked.vstring,
    produceIntermediateHitDoublets = cms.bool(True),
    produceSeedingHitSets = cms.bool(False),
    seedingLayers = cms.InputTag("displacedGeneralStepSeedLayers"),
    trackingRegions = cms.InputTag("displacedGeneralStepTrackingRegions"),
    trackingRegionsSeedingLayers = cms.InputTag("")
)


process.displacedGeneralStepHitTriplets = cms.EDProducer("MultiHitFromChi2EDProducer",
    ClusterShapeHitFilterName = cms.string('ClusterShapeHitFilter'),
    TTRHBuilder = cms.string('WithTrackAngle'),
    chi2VsPtCut = cms.bool(True),
    chi2_cuts = cms.vdouble(3, 4, 5, 5),
    detIdsToDebug = cms.vint32(0, 0, 0),
    doublets = cms.InputTag("displacedGeneralStepHitDoublets"),
    extraHitRPhitolerance = cms.double(0),
    extraHitRZtolerance = cms.double(0),
    extraPhiKDBox = cms.double(0.01),
    extraRKDBox = cms.double(0.2),
    extraZKDBox = cms.double(0.2),
    fnSigmaRZ = cms.double(2),
    maxChi2 = cms.double(5),
    maxElement = cms.uint32(1000000),
    mightGet = cms.optional.untracked.vstring,
    phiPreFiltering = cms.double(0.3),
    pt_interv = cms.vdouble(0.4, 0.7, 1, 2),
    refitHits = cms.bool(True),
    useFixedPreFiltering = cms.bool(False)
)


process.displacedGeneralStepSeedLayers = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(

    ),
    FPix = cms.PSet(

    ),
    MTEC = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTight')
        ),
        maxRing = cms.int32(7),
        minRing = cms.int32(6),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
        skipClusters = cms.InputTag("displacedGeneralStepClusters"),
        useRingSlector = cms.bool(True)
    ),
    MTIB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTight')
        ),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
        skipClusters = cms.InputTag("displacedGeneralStepClusters")
    ),
    MTID = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTight')
        ),
        maxRing = cms.int32(3),
        minRing = cms.int32(3),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
        skipClusters = cms.InputTag("displacedGeneralStepClusters"),
        useRingSlector = cms.bool(True)
    ),
    MTOB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTight')
        ),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
        skipClusters = cms.InputTag("displacedGeneralStepClusters")
    ),
    TEC = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTight')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        maxRing = cms.int32(5),
        minRing = cms.int32(5),
        skipClusters = cms.InputTag("displacedGeneralStepClusters"),
        useRingSlector = cms.bool(True)
    ),
    TIB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTight')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        skipClusters = cms.InputTag("displacedGeneralStepClusters")
    ),
    TID = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTight')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        maxRing = cms.int32(2),
        minRing = cms.int32(1),
        skipClusters = cms.InputTag("displacedGeneralStepClusters"),
        useRingSlector = cms.bool(True)
    ),
    TOB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTight')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        skipClusters = cms.InputTag("displacedGeneralStepClusters")
    ),
    layerList = cms.vstring(
        'TIB1+TIB2+MTIB3',
        'TIB1+TIB2+MTIB4',
        'TIB1+MTIB3+MTIB4',
        'TIB2+MTIB3+MTIB4',
        'TOB1+TOB2+MTOB3',
        'TOB2+MTOB3+MTOB4',
        'MTOB3+MTOB4+MTOB5',
        'MTOB4+MTOB5+MTOB6',
        'MTIB4+TOB1+TOB2',
        'MTIB4+TOB2+MTOB3',
        'MTIB3+TOB1+TOB2',
        'MTID1_pos+TOB1+TOB2',
        'MTID1_neg+TOB1+TOB2',
        'MTID1_pos+TOB1+TOB2',
        'MTID1_neg+TOB1+TOB2',
        'MTID2_pos+TOB1+TOB2',
        'MTID2_neg+TOB1+TOB2',
        'MTID3_pos+TOB1+TOB2',
        'MTID3_neg+TOB1+TOB2',
        'TOB1+TOB2+MTEC1_pos',
        'TOB1+TOB2+MTEC1_neg',
        'TID1+TID2+TEC1_pos',
        'TID1+TID2+TEC1_neg',
        'TID2+MTID3+TEC1_pos',
        'TID2+MTID3+TEC1_neg',
        'MTID3+TEC1_pos+MTEC2_pos',
        'MTID3+TEC1_neg+MTEC2_neg'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.displacedGeneralStepSeeds = cms.EDProducer("SeedCreatorFromRegionConsecutiveHitsTripletOnlyEDProducer",
    MinOneOverPtError = cms.double(1),
    OriginTransverseErrorMultiplier = cms.double(1),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('CombinedSeedComparitor'),
        comparitors = cms.VPSet(
            cms.PSet(
                ClusterShapeCacheSrc = cms.InputTag("siPixelClusterShapeCache"),
                ClusterShapeHitFilterName = cms.string('displacedGeneralStepClusterShapeHitFilter'),
                ComponentName = cms.string('PixelClusterShapeSeedComparitor'),
                FilterAtHelixStage = cms.bool(True),
                FilterPixelHits = cms.bool(False),
                FilterStripHits = cms.bool(True)
            ),
            cms.PSet(
                ComponentName = cms.string('StripSubClusterShapeSeedFilter'),
                FilterAtHelixStage = cms.bool(False),
                label = cms.untracked.string('Seeds'),
                layerMask = cms.PSet(

                ),
                maxNSat = cms.uint32(3),
                maxTrimmedSizeDiffNeg = cms.double(1.0),
                maxTrimmedSizeDiffPos = cms.double(0.7),
                seedCutMIPs = cms.double(0.35),
                seedCutSN = cms.double(7.0),
                subclusterCutMIPs = cms.double(0.45),
                subclusterCutSN = cms.double(12.0),
                subclusterWindow = cms.double(0.7),
                trimMaxADC = cms.double(30.0),
                trimMaxFracNeigh = cms.double(0.25),
                trimMaxFracTotal = cms.double(0.15)
            )
        ),
        mode = cms.string('and')
    ),
    SeedMomentumForBOFF = cms.double(5),
    TTRHBuilder = cms.string('WithTrackAngle'),
    forceKinematicWithRegionDirection = cms.bool(False),
    magneticField = cms.string('ParabolicMf'),
    mightGet = cms.optional.untracked.vstring,
    propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    seedingHitSets = cms.InputTag("displacedGeneralStepHitTriplets")
)


process.displacedGeneralStepTrackCandidates = cms.EDProducer("CkfTrackCandidateMaker",
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    RedundantSeedCleaner = cms.string('CachingSeedCleanerBySharedInput'),
    TrajectoryBuilderPSet = cms.PSet(
        refToPSet_ = cms.string('displacedGeneralStepTrajectoryBuilder')
    ),
    TrajectoryCleaner = cms.string('displacedGeneralStepTrajectoryCleanerBySharedHits'),
    TransientInitialStateEstimatorParameters = cms.PSet(
        numberMeasurementsForFit = cms.int32(4),
        propagatorAlongTISE = cms.string('PropagatorWithMaterial'),
        propagatorOppositeTISE = cms.string('PropagatorWithMaterialOpposite')
    ),
    cleanTrajectoryAfterInOut = cms.bool(True),
    clustersToSkip = cms.InputTag("displacedGeneralStepClusters"),
    doSeedingRegionRebuilding = cms.bool(True),
    maxNSeeds = cms.uint32(500000),
    maxSeedsBeforeCleaning = cms.uint32(5000),
    mightGet = cms.optional.untracked.vstring,
    numHitsForSeedCleaner = cms.int32(50),
    onlyPixelHitsForSeedCleaner = cms.bool(False),
    phase2clustersToSkip = cms.InputTag(""),
    reverseTrajectories = cms.bool(False),
    src = cms.InputTag("displacedGeneralStepSeeds"),
    useHitsSplitting = cms.bool(True)
)


process.displacedGeneralStepTrackingRegions = cms.EDProducer("GlobalTrackingRegionEDProducer",
    RegionPSet = cms.PSet(
        originHalfLength = cms.double(55),
        originRadius = cms.double(10),
        originXPos = cms.double(0),
        originYPos = cms.double(0),
        originZPos = cms.double(0),
        precise = cms.bool(True),
        ptMin = cms.double(1),
        useMultipleScattering = cms.bool(True)
    ),
    mightGet = cms.optional.untracked.vstring
)


process.displacedGeneralStepTracks = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('displacedGeneralStep'),
    Fitter = cms.string('generalDisplacedFlexibleKFFittingSmoother'),
    GeometricInnerState = cms.bool(False),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    Propagator = cms.string('RungeKuttaTrackerPropagator'),
    SimpleMagneticField = cms.string(''),
    TTRHBuilder = cms.string('WithAngleAndTemplate'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("displacedGeneralStepTrackCandidates"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(False)
)


process.duplicateTrackCandidates = cms.EDProducer("DuplicateTrackMerger",
    GBRForestFileName = cms.string(''),
    chi2EstimatorName = cms.string('duplicateTrackCandidatesChi2Est'),
    forestLabel = cms.string('MVADuplicate'),
    maxDCA = cms.double(30),
    maxDLambda = cms.double(0.3),
    maxDPhi = cms.double(0.3),
    maxDQoP = cms.double(0.25),
    maxDdsz = cms.double(10),
    maxDdxy = cms.double(10),
    mightGet = cms.optional.untracked.vstring,
    minBDTG = cms.double(-0.1),
    minDeltaR3d = cms.double(-4),
    minP = cms.double(0.4),
    minpT = cms.double(0.2),
    overlapCheckMaxHits = cms.uint32(4),
    overlapCheckMaxMissingLayers = cms.uint32(1),
    overlapCheckMinCosT = cms.double(0.99),
    propagatorName = cms.string('PropagatorWithMaterial'),
    source = cms.InputTag("preDuplicateMergingGeneralTracks"),
    ttrhBuilderName = cms.string('WithAngleAndTemplate'),
    useInnermostState = cms.bool(True)
)


process.duplicateTrackClassifier = cms.EDProducer("TrackCutClassifier",
    beamspot = cms.InputTag("offlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mightGet = cms.optional.untracked.vstring,
    mva = cms.PSet(
        dr_par = cms.PSet(
            d0err = cms.vdouble(0.003, 0.003, 0.003),
            d0err_par = cms.vdouble(0.001, 0.001, 0.001),
            drWPVerr_par = cms.vdouble(3.4028234663852886e+38, 3.4028234663852886e+38, 3.4028234663852886e+38),
            dr_exp = cms.vint32(2147483647, 2147483647, 2147483647),
            dr_par1 = cms.vdouble(3.4028234663852886e+38, 3.4028234663852886e+38, 3.4028234663852886e+38),
            dr_par2 = cms.vdouble(3.4028234663852886e+38, 3.4028234663852886e+38, 3.4028234663852886e+38)
        ),
        dz_par = cms.PSet(
            dzWPVerr_par = cms.vdouble(3.4028234663852886e+38, 3.4028234663852886e+38, 3.4028234663852886e+38),
            dz_exp = cms.vint32(2147483647, 2147483647, 2147483647),
            dz_par1 = cms.vdouble(3.4028234663852886e+38, 3.4028234663852886e+38, 3.4028234663852886e+38),
            dz_par2 = cms.vdouble(3.4028234663852886e+38, 3.4028234663852886e+38, 3.4028234663852886e+38)
        ),
        isHLT = cms.bool(False),
        maxChi2 = cms.vdouble(9999.0, 9999.0, 9999.0),
        maxChi2n = cms.vdouble(10.0, 1.0, 0.4),
        maxDr = cms.vdouble(3.4028234663852886e+38, 3.4028234663852886e+38, 3.4028234663852886e+38),
        maxDz = cms.vdouble(3.4028234663852886e+38, 3.4028234663852886e+38, 3.4028234663852886e+38),
        maxDzWrtBS = cms.vdouble(3.4028234663852886e+38, 24, 15),
        maxLostLayers = cms.vint32(99, 99, 99),
        maxRelPtErr = cms.vdouble(3.4028234663852886e+38, 3.4028234663852886e+38, 3.4028234663852886e+38),
        min3DLayers = cms.vint32(0, 0, 0),
        minHits = cms.vint32(0, 0, 1),
        minHits4pass = cms.vint32(2147483647, 2147483647, 2147483647),
        minLayers = cms.vint32(0, 0, 0),
        minNVtxTrk = cms.int32(2),
        minNdof = cms.vdouble(-1, -1, -1),
        minPixelHits = cms.vint32(0, 0, 0)
    ),
    qualityCuts = cms.vdouble(-0.7, 0.1, 0.7),
    src = cms.InputTag("mergedDuplicateTracks"),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.earlyGeneralTracks = cms.EDProducer("TrackCollectionMerger",
    allowFirstHitShare = cms.bool(True),
    copyExtras = cms.untracked.bool(True),
    copyTrajectories = cms.untracked.bool(False),
    enableMerging = cms.bool(True),
    foundHitBonus = cms.double(10),
    inputClassifiers = cms.vstring(
        'initialStep',
        'highPtTripletStep',
        'jetCoreRegionalStep',
        'lowPtQuadStep',
        'lowPtTripletStep',
        'detachedQuadStep',
        'detachedTripletStep',
        'pixelPairStep',
        'mixedTripletStep',
        'pixelLessStep',
        'tobTecStep'
    ),
    lostHitPenalty = cms.double(5),
    mightGet = cms.optional.untracked.vstring,
    minQuality = cms.string('loose'),
    minShareHits = cms.uint32(2),
    shareFrac = cms.double(0.19),
    trackAlgoPriorityOrder = cms.string('trackAlgoPriorityOrder'),
    trackProducers = cms.VInputTag(
        "initialStepTracks", "highPtTripletStepTracks", "jetCoreRegionalStepTracks", "lowPtQuadStepTracks", "lowPtTripletStepTracks",
        "detachedQuadStepTracks", "detachedTripletStepTracks", "pixelPairStepTracks", "mixedTripletStepTracks", "pixelLessStepTracks",
        "tobTecStepTracks"
    )
)


process.earlyMuons = cms.EDProducer("MuonIdProducer",
    CaloExtractorPSet = cms.PSet(
        CenterConeOnCalIntersection = cms.bool(False),
        ComponentName = cms.string('CaloExtractorByAssociator'),
        DR_Max = cms.double(0.5),
        DR_Veto_E = cms.double(0.07),
        DR_Veto_H = cms.double(0.1),
        DR_Veto_HO = cms.double(0.1),
        DepositInstanceLabels = cms.vstring(
            'ecal',
            'hcal',
            'ho'
        ),
        DepositLabel = cms.untracked.string('Cal'),
        NoiseTow_EB = cms.double(0.04),
        NoiseTow_EE = cms.double(0.15),
        Noise_EB = cms.double(0.025),
        Noise_EE = cms.double(0.1),
        Noise_HB = cms.double(0.2),
        Noise_HE = cms.double(0.2),
        Noise_HO = cms.double(0.2),
        PrintTimeReport = cms.untracked.bool(False),
        PropagatorName = cms.string('SteppingHelixPropagatorAny'),
        ServiceParameters = cms.PSet(
            Propagators = cms.untracked.vstring('SteppingHelixPropagatorAny'),
            RPCLayers = cms.bool(False),
            UseMuonNavigation = cms.untracked.bool(False)
        ),
        Threshold_E = cms.double(0.2),
        Threshold_H = cms.double(0.5),
        Threshold_HO = cms.double(0.5),
        TrackAssociatorParameters = cms.PSet(
            CSCSegmentCollectionLabel = cms.InputTag("cscSegments"),
            CaloTowerCollectionLabel = cms.InputTag("towerMaker"),
            DTRecSegment4DCollectionLabel = cms.InputTag("dt4DSegments"),
            EBRecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
            EERecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
            GEMSegmentCollectionLabel = cms.InputTag("gemSegments"),
            HBHERecHitCollectionLabel = cms.InputTag("hbhereco"),
            HORecHitCollectionLabel = cms.InputTag("horeco"),
            ME0SegmentCollectionLabel = cms.InputTag("me0Segments"),
            accountForTrajectoryChangeCalo = cms.bool(False),
            dREcal = cms.double(1.0),
            dREcalPreselection = cms.double(1.0),
            dRHcal = cms.double(1.0),
            dRHcalPreselection = cms.double(1.0),
            dRMuon = cms.double(9999.0),
            dRMuonPreselection = cms.double(0.2),
            dRPreshowerPreselection = cms.double(0.2),
            muonMaxDistanceSigmaX = cms.double(0.0),
            muonMaxDistanceSigmaY = cms.double(0.0),
            muonMaxDistanceX = cms.double(5.0),
            muonMaxDistanceY = cms.double(5.0),
            propagateAllDirections = cms.bool(True),
            trajectoryUncertaintyTolerance = cms.double(-1.0),
            truthMatch = cms.bool(False),
            useCalo = cms.bool(True),
            useEcal = cms.bool(False),
            useGEM = cms.bool(False),
            useHO = cms.bool(False),
            useHcal = cms.bool(False),
            useME0 = cms.bool(False),
            useMuon = cms.bool(False),
            usePreshower = cms.bool(False)
        ),
        UseRecHitsFlag = cms.bool(False)
    ),
    JetExtractorPSet = cms.PSet(
        ComponentName = cms.string('JetExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(0.1),
        ExcludeMuonVeto = cms.bool(True),
        JetCollectionLabel = cms.InputTag("ak4CaloJets"),
        PrintTimeReport = cms.untracked.bool(False),
        PropagatorName = cms.string('SteppingHelixPropagatorAny'),
        ServiceParameters = cms.PSet(
            Propagators = cms.untracked.vstring('SteppingHelixPropagatorAny'),
            RPCLayers = cms.bool(False),
            UseMuonNavigation = cms.untracked.bool(False)
        ),
        Threshold = cms.double(5.0),
        TrackAssociatorParameters = cms.PSet(
            CSCSegmentCollectionLabel = cms.InputTag("cscSegments"),
            CaloTowerCollectionLabel = cms.InputTag("towerMaker"),
            DTRecSegment4DCollectionLabel = cms.InputTag("dt4DSegments"),
            EBRecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
            EERecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
            GEMSegmentCollectionLabel = cms.InputTag("gemSegments"),
            HBHERecHitCollectionLabel = cms.InputTag("hbhereco"),
            HORecHitCollectionLabel = cms.InputTag("horeco"),
            ME0SegmentCollectionLabel = cms.InputTag("me0Segments"),
            accountForTrajectoryChangeCalo = cms.bool(False),
            dREcal = cms.double(0.5),
            dREcalPreselection = cms.double(0.5),
            dRHcal = cms.double(0.5),
            dRHcalPreselection = cms.double(0.5),
            dRMuon = cms.double(9999.0),
            dRMuonPreselection = cms.double(0.2),
            dRPreshowerPreselection = cms.double(0.2),
            muonMaxDistanceSigmaX = cms.double(0.0),
            muonMaxDistanceSigmaY = cms.double(0.0),
            muonMaxDistanceX = cms.double(5.0),
            muonMaxDistanceY = cms.double(5.0),
            propagateAllDirections = cms.bool(True),
            trajectoryUncertaintyTolerance = cms.double(-1.0),
            truthMatch = cms.bool(False),
            useCalo = cms.bool(True),
            useEcal = cms.bool(False),
            useGEM = cms.bool(False),
            useHO = cms.bool(False),
            useHcal = cms.bool(False),
            useME0 = cms.bool(False),
            useMuon = cms.bool(False),
            usePreshower = cms.bool(False)
        )
    ),
    MuonCaloCompatibility = cms.PSet(
        MuonTemplateFileName = cms.FileInPath('RecoMuon/MuonIdentification/data/MuID_templates_muons_lowPt_3_1_norm.root'),
        PionTemplateFileName = cms.FileInPath('RecoMuon/MuonIdentification/data/MuID_templates_pions_lowPt_3_1_norm.root'),
        allSiPMHO = cms.bool(False),
        delta_eta = cms.double(0.02),
        delta_phi = cms.double(0.02)
    ),
    ShowerDigiFillerParameters = cms.PSet(
        cscDigiCollectionLabel = cms.InputTag("muonCSCDigis","MuonCSCStripDigi"),
        digiMaxDistanceX = cms.double(25.0),
        dtDigiCollectionLabel = cms.InputTag("muonDTDigis")
    ),
    TimingFillerParameters = cms.PSet(
        CSCTimingParameters = cms.PSet(
            CSCStripError = cms.double(7.0),
            CSCStripTimeOffset = cms.double(0.0),
            CSCWireError = cms.double(8.6),
            CSCWireTimeOffset = cms.double(0.0),
            PruneCut = cms.double(9.0),
            ServiceParameters = cms.PSet(
                Propagators = cms.untracked.vstring(
                    'SteppingHelixPropagatorAny',
                    'PropagatorWithMaterial',
                    'PropagatorWithMaterialOpposite'
                ),
                RPCLayers = cms.bool(True)
            ),
            UseStripTime = cms.bool(True),
            UseWireTime = cms.bool(True),
            debug = cms.bool(False)
        ),
        DTTimingParameters = cms.PSet(
            DTTimeOffset = cms.double(0.0),
            DoWireCorr = cms.bool(True),
            DropTheta = cms.bool(True),
            HitError = cms.double(2.8),
            HitsMin = cms.int32(3),
            PruneCut = cms.double(5.0),
            RequireBothProjections = cms.bool(False),
            ServiceParameters = cms.PSet(
                Propagators = cms.untracked.vstring(
                    'SteppingHelixPropagatorAny',
                    'PropagatorWithMaterial',
                    'PropagatorWithMaterialOpposite'
                ),
                RPCLayers = cms.bool(True)
            ),
            UseSegmentT0 = cms.bool(False),
            debug = cms.bool(False)
        ),
        EcalEnergyCut = cms.double(0.4),
        ErrorEB = cms.double(2.085),
        ErrorEE = cms.double(6.95),
        MatchParameters = cms.PSet(
            CSCsegments = cms.InputTag("cscSegments"),
            DTradius = cms.double(0.01),
            DTsegments = cms.InputTag("dt4DSegments"),
            RPChits = cms.InputTag("rpcRecHits"),
            TightMatchCSC = cms.bool(True),
            TightMatchDT = cms.bool(False)
        ),
        UseCSC = cms.bool(True),
        UseDT = cms.bool(True),
        UseECAL = cms.bool(False)
    ),
    TrackAssociatorParameters = cms.PSet(
        CSCSegmentCollectionLabel = cms.InputTag("cscSegments"),
        CaloTowerCollectionLabel = cms.InputTag("towerMaker"),
        DTRecSegment4DCollectionLabel = cms.InputTag("dt4DSegments"),
        EBRecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
        EERecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
        GEMSegmentCollectionLabel = cms.InputTag("gemSegments"),
        HBHERecHitCollectionLabel = cms.InputTag("hbhereco"),
        HORecHitCollectionLabel = cms.InputTag("horeco"),
        ME0SegmentCollectionLabel = cms.InputTag("me0Segments"),
        accountForTrajectoryChangeCalo = cms.bool(False),
        dREcal = cms.double(9999.0),
        dREcalPreselection = cms.double(0.05),
        dRHcal = cms.double(9999.0),
        dRHcalPreselection = cms.double(0.2),
        dRMuon = cms.double(9999.0),
        dRMuonPreselection = cms.double(0.2),
        dRPreshowerPreselection = cms.double(0.2),
        muonMaxDistanceSigmaX = cms.double(0.0),
        muonMaxDistanceSigmaY = cms.double(0.0),
        muonMaxDistanceX = cms.double(5.0),
        muonMaxDistanceY = cms.double(5.0),
        propagateAllDirections = cms.bool(True),
        trajectoryUncertaintyTolerance = cms.double(-1.0),
        truthMatch = cms.bool(False),
        useCalo = cms.bool(False),
        useEcal = cms.bool(False),
        useGEM = cms.bool(True),
        useHO = cms.bool(False),
        useHcal = cms.bool(False),
        useME0 = cms.bool(False),
        useMuon = cms.bool(True),
        usePreshower = cms.bool(False)
    ),
    TrackExtractorPSet = cms.PSet(
        BeamSpotLabel = cms.InputTag("offlineBeamSpot"),
        BeamlineOption = cms.string('BeamSpotFromEvent'),
        Chi2Ndof_Max = cms.double(1e+64),
        Chi2Prob_Min = cms.double(-1.0),
        ComponentName = cms.string('TrackExtractor'),
        DR_Max = cms.double(0.5),
        DR_Veto = cms.double(0.01),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(0.1),
        Diff_z = cms.double(0.2),
        NHits_Min = cms.uint32(0),
        Pt_Min = cms.double(-1.0),
        inputTrackCollection = cms.InputTag("generalTracks")
    ),
    TrackerKinkFinderParameters = cms.PSet(
        DoPredictionsOnly = cms.bool(False),
        Fitter = cms.string('KFFitterForRefitInsideOut'),
        MTDRecHitBuilder = cms.string('MTDRecHitBuilder'),
        MuonRecHitBuilder = cms.string('MuonRecHitBuilder'),
        Propagator = cms.string('SmartPropagatorAnyRKOpposite'),
        RefitDirection = cms.string('alongMomentum'),
        RefitRPCHits = cms.bool(True),
        Smoother = cms.string('KFSmootherForRefitInsideOut'),
        TrackerRecHitBuilder = cms.string('WithAngleAndTemplate'),
        diagonalOnly = cms.bool(False),
        usePosition = cms.bool(True)
    ),
    addExtraSoftMuons = cms.bool(False),
    arbitrateTrackerMuons = cms.bool(True),
    arbitrationCleanerOptions = cms.PSet(
        ClusterDPhi = cms.double(0.6),
        ClusterDTheta = cms.double(0.02),
        Clustering = cms.bool(True),
        ME1a = cms.bool(True),
        Overlap = cms.bool(True),
        OverlapDPhi = cms.double(0.0786),
        OverlapDTheta = cms.double(0.02)
    ),
    debugWithTruthMatching = cms.bool(False),
    ecalDepositName = cms.string('ecal'),
    fillCaloCompatibility = cms.bool(False),
    fillEnergy = cms.bool(False),
    fillGlobalTrackQuality = cms.bool(False),
    fillGlobalTrackRefits = cms.bool(False),
    fillIsolation = cms.bool(False),
    fillMatching = cms.bool(True),
    fillShowerDigis = cms.bool(True),
    fillTrackerKink = cms.bool(False),
    globalTrackQualityInputTag = cms.InputTag("glbTrackQual"),
    hcalDepositName = cms.string('hcal'),
    hoDepositName = cms.string('ho'),
    inputCollectionLabels = cms.VInputTag("earlyGeneralTracks", "standAloneMuons:UpdatedAtVtx"),
    inputCollectionTypes = cms.vstring(
        'inner tracks',
        'outer tracks'
    ),
    jetDepositName = cms.string('jets'),
    maxAbsDx = cms.double(3.0),
    maxAbsDy = cms.double(9999.0),
    maxAbsEta = cms.double(3.0),
    maxAbsPullX = cms.double(3.0),
    maxAbsPullY = cms.double(9999.0),
    minCaloCompatibility = cms.double(0.6),
    minNumberOfMatches = cms.int32(1),
    minP = cms.double(3.0),
    minPCaloMuon = cms.double(3.0),
    minPt = cms.double(2.0),
    ptThresholdToFillCandidateP4WithGlobalFit = cms.double(200.0),
    pvInputTag = cms.InputTag("offlinePrimaryVertices"),
    runArbitrationCleaner = cms.bool(True),
    selectHighPurity = cms.bool(False),
    sigmaThresholdToFillCandidateP4WithGlobalFit = cms.double(2.0),
    storeCrossedHcalRecHits = cms.bool(True),
    trackDepositName = cms.string('tracker'),
    writeIsoDeposits = cms.bool(True)
)


process.firstStepPrimaryVertices = cms.EDProducer("RecoChargedRefCandidatePrimaryVertexSorter",
    assignment = cms.PSet(
        DzCutForChargedFromPUVtxs = cms.double(0.2),
        EtaMinUseDz = cms.double(-1),
        NumOfPUVtxsForCharged = cms.uint32(0),
        OnlyUseFirstDz = cms.bool(False),
        PtMaxCharged = cms.double(-1),
        maxDistanceToJetAxis = cms.double(0.07),
        maxDtSigForPrimaryAssignment = cms.double(3),
        maxDxyForJetAxisAssigment = cms.double(0.1),
        maxDxyForNotReconstructedPrimary = cms.double(0.01),
        maxDxySigForNotReconstructedPrimary = cms.double(2),
        maxDzErrorForPrimaryAssignment = cms.double(0.05),
        maxDzForJetAxisAssigment = cms.double(0.1),
        maxDzForPrimaryAssignment = cms.double(0.1),
        maxDzSigForPrimaryAssignment = cms.double(5),
        maxJetDeltaR = cms.double(0.5),
        minJetPt = cms.double(25),
        preferHighRanked = cms.bool(False),
        useTiming = cms.bool(False),
        useVertexFit = cms.bool(True)
    ),
    jets = cms.InputTag("ak4CaloJetsForTrk"),
    mightGet = cms.optional.untracked.vstring,
    particles = cms.InputTag("initialStepTrackRefsForJets"),
    produceAssociationToOriginalVertices = cms.bool(False),
    produceNoPileUpCollection = cms.bool(False),
    producePileUpCollection = cms.bool(False),
    produceSortedVertices = cms.bool(True),
    qualityForPrimary = cms.int32(3),
    sorting = cms.PSet(

    ),
    trackTimeResoTag = cms.InputTag(""),
    trackTimeTag = cms.InputTag(""),
    usePVMET = cms.bool(True),
    vertices = cms.InputTag("firstStepPrimaryVerticesUnsorted")
)


process.firstStepPrimaryVerticesBeforeMixing = cms.EDProducer("PrimaryVertexProducer",
    TkClusParameters = cms.PSet(
        TkDAClusParameters = cms.PSet(
            Tmin = cms.double(2.0),
            Tpurge = cms.double(2.0),
            Tstop = cms.double(0.5),
            block_size = cms.uint32(10000),
            convergence_mode = cms.int32(0),
            coolingFactor = cms.double(0.6),
            d0CutOff = cms.double(3.0),
            delta_highT = cms.double(0.01),
            delta_lowT = cms.double(0.001),
            dzCutOff = cms.double(3.0),
            overlap_frac = cms.double(0.0),
            runInBlocks = cms.bool(False),
            uniquetrkminp = cms.double(0.0),
            uniquetrkweight = cms.double(0.8),
            vertexSize = cms.double(0.006),
            zmerge = cms.double(0.01),
            zrange = cms.double(4.0)
        ),
        algorithm = cms.string('DA_vect')
    ),
    TkFilterParameters = cms.PSet(
        algorithm = cms.string('filter'),
        maxD0Error = cms.double(1.0),
        maxD0Significance = cms.double(4.0),
        maxDzError = cms.double(1.0),
        maxEta = cms.double(2.4),
        maxNormalizedChi2 = cms.double(10.0),
        minPixelLayersWithHits = cms.int32(2),
        minPt = cms.double(0.0),
        minSiliconLayersWithHits = cms.int32(5),
        trackQuality = cms.string('any')
    ),
    TrackLabel = cms.InputTag("initialStepTracks"),
    beamSpotLabel = cms.InputTag("offlineBeamSpot"),
    isRecoveryIteration = cms.bool(False),
    recoveryVtxCollection = cms.InputTag(""),
    verbose = cms.untracked.bool(False),
    vertexCollections = cms.VPSet(cms.PSet(
        algorithm = cms.string('AdaptiveVertexFitter'),
        chi2cutoff = cms.double(2.5),
        label = cms.string(''),
        maxDistanceToBeam = cms.double(1.0),
        minNdof = cms.double(0.0),
        useBeamConstraint = cms.bool(False)
    ))
)


process.firstStepPrimaryVerticesPreSplitting = cms.EDProducer("PrimaryVertexProducer",
    TkClusParameters = cms.PSet(
        TkDAClusParameters = cms.PSet(
            Tmin = cms.double(2.0),
            Tpurge = cms.double(2.0),
            Tstop = cms.double(0.5),
            block_size = cms.uint32(10000),
            convergence_mode = cms.int32(0),
            coolingFactor = cms.double(0.6),
            d0CutOff = cms.double(3.0),
            delta_highT = cms.double(0.01),
            delta_lowT = cms.double(0.001),
            dzCutOff = cms.double(3.0),
            overlap_frac = cms.double(0.0),
            runInBlocks = cms.bool(False),
            uniquetrkminp = cms.double(0.0),
            uniquetrkweight = cms.double(0.8),
            vertexSize = cms.double(0.006),
            zmerge = cms.double(0.01),
            zrange = cms.double(4.0)
        ),
        algorithm = cms.string('DA_vect')
    ),
    TkFilterParameters = cms.PSet(
        algorithm = cms.string('filter'),
        maxD0Error = cms.double(1.0),
        maxD0Significance = cms.double(4.0),
        maxDzError = cms.double(1.0),
        maxEta = cms.double(2.4),
        maxNormalizedChi2 = cms.double(10.0),
        minPixelLayersWithHits = cms.int32(2),
        minPt = cms.double(0.0),
        minSiliconLayersWithHits = cms.int32(5),
        trackQuality = cms.string('any')
    ),
    TrackLabel = cms.InputTag("initialStepTracksPreSplitting"),
    beamSpotLabel = cms.InputTag("offlineBeamSpot"),
    isRecoveryIteration = cms.bool(False),
    recoveryVtxCollection = cms.InputTag(""),
    verbose = cms.untracked.bool(False),
    vertexCollections = cms.VPSet(cms.PSet(
        algorithm = cms.string('AdaptiveVertexFitter'),
        chi2cutoff = cms.double(2.5),
        label = cms.string(''),
        maxDistanceToBeam = cms.double(1.0),
        minNdof = cms.double(0.0),
        useBeamConstraint = cms.bool(False)
    ))
)


process.firstStepPrimaryVerticesUnsorted = cms.EDProducer("PrimaryVertexProducer",
    TkClusParameters = cms.PSet(
        TkDAClusParameters = cms.PSet(
            Tmin = cms.double(2.0),
            Tpurge = cms.double(2.0),
            Tstop = cms.double(0.5),
            block_size = cms.uint32(10000),
            convergence_mode = cms.int32(0),
            coolingFactor = cms.double(0.6),
            d0CutOff = cms.double(3.0),
            delta_highT = cms.double(0.01),
            delta_lowT = cms.double(0.001),
            dzCutOff = cms.double(3.0),
            overlap_frac = cms.double(0.0),
            runInBlocks = cms.bool(False),
            uniquetrkminp = cms.double(0.0),
            uniquetrkweight = cms.double(0.8),
            vertexSize = cms.double(0.006),
            zmerge = cms.double(0.01),
            zrange = cms.double(4.0)
        ),
        algorithm = cms.string('DA_vect')
    ),
    TkFilterParameters = cms.PSet(
        algorithm = cms.string('filter'),
        maxD0Error = cms.double(1.0),
        maxD0Significance = cms.double(4.0),
        maxDzError = cms.double(1.0),
        maxEta = cms.double(2.4),
        maxNormalizedChi2 = cms.double(10.0),
        minPixelLayersWithHits = cms.int32(2),
        minPt = cms.double(0.0),
        minSiliconLayersWithHits = cms.int32(5),
        trackQuality = cms.string('any')
    ),
    TrackLabel = cms.InputTag("initialStepTracks"),
    beamSpotLabel = cms.InputTag("offlineBeamSpot"),
    isRecoveryIteration = cms.bool(False),
    recoveryVtxCollection = cms.InputTag(""),
    verbose = cms.untracked.bool(False),
    vertexCollections = cms.VPSet(cms.PSet(
        algorithm = cms.string('AdaptiveVertexFitter'),
        chi2cutoff = cms.double(2.5),
        label = cms.string(''),
        maxDistanceToBeam = cms.double(1.0),
        minNdof = cms.double(0.0),
        useBeamConstraint = cms.bool(False)
    ))
)


process.generalTracks = cms.EDProducer("DuplicateListMerger",
    candidateComponents = cms.InputTag("duplicateTrackCandidates","candidateMap"),
    candidateSource = cms.InputTag("duplicateTrackCandidates","candidates"),
    copyExtras = cms.untracked.bool(True),
    copyTrajectories = cms.untracked.bool(False),
    diffHitsCut = cms.int32(5),
    mergedMVAVals = cms.InputTag("duplicateTrackClassifier","MVAValues"),
    mergedSource = cms.InputTag("mergedDuplicateTracks"),
    mightGet = cms.optional.untracked.vstring,
    originalMVAVals = cms.InputTag("preDuplicateMergingGeneralTracks","MVAValues"),
    originalSource = cms.InputTag("preDuplicateMergingGeneralTracks"),
    trackAlgoPriorityOrder = cms.string('trackAlgoPriorityOrder')
)


process.globalCombinedSeeds = cms.EDProducer("SeedCombiner",
    seedCollections = cms.VInputTag(cms.InputTag("globalSeedsFromTripletsWithVertices"), cms.InputTag("globalSeedsFromPairsWithVertices"))
)


process.globalMixedSeeds = cms.EDProducer("SeedGeneratorFromRegionHitsEDProducer",
    ClusterCheckPSet = cms.PSet(
        ClusterCollectionLabel = cms.InputTag("siStripClusters"),
        MaxNumberOfCosmicClusters = cms.uint32(400000),
        MaxNumberOfPixelClusters = cms.uint32(40000),
        PixelClusterCollectionLabel = cms.InputTag("siPixelClusters"),
        cut = cms.string('strip < 400000 && pixel < 40000 && (strip < 50000 + 10*pixel) && (pixel < 5000 + 0.1*strip)'),
        doClusterCheck = cms.bool(False)
    ),
    OrderedHitsFactoryPSet = cms.PSet(
        ComponentName = cms.string('StandardHitPairGenerator'),
        SeedingLayers = cms.InputTag("MixedLayerPairs"),
        maxElement = cms.uint32(1000000)
    ),
    RegionFactoryPSet = cms.PSet(
        ComponentName = cms.string('GlobalRegionProducerFromBeamSpot'),
        RegionPSet = cms.PSet(
            beamSpot = cms.InputTag("offlineBeamSpot"),
            originHalfLength = cms.double(21.2),
            originRadius = cms.double(0.2),
            precise = cms.bool(True),
            ptMin = cms.double(0.9),
            useMultipleScattering = cms.bool(False)
        )
    ),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('none')
    ),
    SeedCreatorPSet = cms.PSet(
        ComponentName = cms.string('SeedFromConsecutiveHitsCreator'),
        MinOneOverPtError = cms.double(1.0),
        OriginTransverseErrorMultiplier = cms.double(1.0),
        SeedMomentumForBOFF = cms.double(5.0),
        TTRHBuilder = cms.string('WithTrackAngle'),
        forceKinematicWithRegionDirection = cms.bool(False),
        magneticField = cms.string('ParabolicMf'),
        propagator = cms.string('PropagatorWithMaterialParabolicMf')
    )
)


process.globalPixelLessSeeds = cms.EDProducer("SeedGeneratorFromRegionHitsEDProducer",
    ClusterCheckPSet = cms.PSet(
        ClusterCollectionLabel = cms.InputTag("siStripClusters"),
        MaxNumberOfCosmicClusters = cms.uint32(5000),
        MaxNumberOfPixelClusters = cms.uint32(40000),
        PixelClusterCollectionLabel = cms.InputTag("siPixelClusters"),
        cut = cms.string('strip < 400000 && pixel < 40000 && (strip < 50000 + 10*pixel) && (pixel < 5000 + 0.1*strip)'),
        doClusterCheck = cms.bool(False)
    ),
    OrderedHitsFactoryPSet = cms.PSet(
        ComponentName = cms.string('StandardHitPairGenerator'),
        SeedingLayers = cms.InputTag("pixelLessLayerPairs4PixelLessTracking"),
        maxElement = cms.uint32(100000)
    ),
    RegionFactoryPSet = cms.PSet(
        ComponentName = cms.string('GlobalRegionProducerFromBeamSpot'),
        RegionPSet = cms.PSet(
            beamSpot = cms.InputTag("offlineBeamSpot"),
            originHalfLength = cms.double(40),
            originRadius = cms.double(0.2),
            precise = cms.bool(True),
            ptMin = cms.double(0.9),
            useMultipleScattering = cms.bool(False)
        )
    ),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('none')
    ),
    SeedCreatorPSet = cms.PSet(
        ComponentName = cms.string('SeedFromConsecutiveHitsCreator'),
        MinOneOverPtError = cms.double(1.0),
        OriginTransverseErrorMultiplier = cms.double(1.0),
        SeedMomentumForBOFF = cms.double(5.0),
        TTRHBuilder = cms.string('WithTrackAngle'),
        forceKinematicWithRegionDirection = cms.bool(False),
        magneticField = cms.string('ParabolicMf'),
        propagator = cms.string('PropagatorWithMaterialParabolicMf')
    )
)


process.globalPixelSeeds = cms.EDProducer("SeedGeneratorFromRegionHitsEDProducer",
    ClusterCheckPSet = cms.PSet(
        ClusterCollectionLabel = cms.InputTag("siStripClusters"),
        MaxNumberOfCosmicClusters = cms.uint32(400000),
        MaxNumberOfPixelClusters = cms.uint32(40000),
        PixelClusterCollectionLabel = cms.InputTag("siPixelClusters"),
        cut = cms.string('strip < 400000 && pixel < 40000 && (strip < 50000 + 10*pixel) && (pixel < 5000 + 0.1*strip)'),
        doClusterCheck = cms.bool(False)
    ),
    OrderedHitsFactoryPSet = cms.PSet(
        ComponentName = cms.string('StandardHitPairGenerator'),
        SeedingLayers = cms.InputTag("PixelLayerPairs")
    ),
    RegionFactoryPSet = cms.PSet(
        ComponentName = cms.string('GlobalRegionProducerFromBeamSpot'),
        RegionPSet = cms.PSet(
            beamSpot = cms.InputTag("offlineBeamSpot"),
            originHalfLength = cms.double(21.2),
            originRadius = cms.double(0.2),
            precise = cms.bool(True),
            ptMin = cms.double(0.9),
            useMultipleScattering = cms.bool(False)
        )
    ),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('none')
    ),
    SeedCreatorPSet = cms.PSet(
        ComponentName = cms.string('SeedFromConsecutiveHitsCreator'),
        MinOneOverPtError = cms.double(1.0),
        OriginTransverseErrorMultiplier = cms.double(1.0),
        SeedMomentumForBOFF = cms.double(5.0),
        TTRHBuilder = cms.string('WithTrackAngle'),
        forceKinematicWithRegionDirection = cms.bool(False),
        magneticField = cms.string('ParabolicMf'),
        propagator = cms.string('PropagatorWithMaterialParabolicMf')
    )
)


process.globalSeedsFromPairsWithVertices = cms.EDProducer("SeedGeneratorFromRegionHitsEDProducer",
    ClusterCheckPSet = cms.PSet(
        ClusterCollectionLabel = cms.InputTag("siStripClusters"),
        MaxNumberOfCosmicClusters = cms.uint32(400000),
        MaxNumberOfPixelClusters = cms.uint32(40000),
        PixelClusterCollectionLabel = cms.InputTag("siPixelClusters"),
        cut = cms.string('strip < 400000 && pixel < 40000 && (strip < 50000 + 10*pixel) && (pixel < 5000 + 0.1*strip)'),
        doClusterCheck = cms.bool(False)
    ),
    OrderedHitsFactoryPSet = cms.PSet(
        ComponentName = cms.string('StandardHitPairGenerator'),
        SeedingLayers = cms.InputTag("MixedLayerPairs"),
        maxElement = cms.uint32(1000000)
    ),
    RegionFactoryPSet = cms.PSet(
        ComponentName = cms.string('GlobalTrackingRegionWithVerticesProducer'),
        RegionPSet = cms.PSet(
            VertexCollection = cms.InputTag("firstStepPrimaryVertices"),
            beamSpot = cms.InputTag("offlineBeamSpot"),
            fixedError = cms.double(0.2),
            halfLengthScaling4BigEvts = cms.bool(False),
            maxNVertices = cms.int32(-1),
            maxPtMin = cms.double(1000),
            minHalfLength = cms.double(0),
            minOriginR = cms.double(0),
            nSigmaZ = cms.double(4),
            originRScaling4BigEvts = cms.bool(False),
            originRadius = cms.double(0.2),
            pixelClustersForScaling = cms.InputTag("siPixelClusters"),
            precise = cms.bool(True),
            ptMin = cms.double(0.9),
            ptMinScaling4BigEvts = cms.bool(False),
            scalingEndNPix = cms.double(1),
            scalingStartNPix = cms.double(0),
            sigmaZVertex = cms.double(3),
            useFakeVertices = cms.bool(False),
            useFixedError = cms.bool(True),
            useFoundVertices = cms.bool(True),
            useMultipleScattering = cms.bool(False)
        )
    ),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('none')
    ),
    SeedCreatorPSet = cms.PSet(
        ComponentName = cms.string('SeedFromConsecutiveHitsCreator'),
        MinOneOverPtError = cms.double(1.0),
        OriginTransverseErrorMultiplier = cms.double(1.0),
        SeedMomentumForBOFF = cms.double(5.0),
        TTRHBuilder = cms.string('WithTrackAngle'),
        forceKinematicWithRegionDirection = cms.bool(False),
        magneticField = cms.string('ParabolicMf'),
        propagator = cms.string('PropagatorWithMaterialParabolicMf')
    )
)


process.globalSeedsFromTriplets = cms.EDProducer("SeedGeneratorFromRegionHitsEDProducer",
    ClusterCheckPSet = cms.PSet(
        ClusterCollectionLabel = cms.InputTag("siStripClusters"),
        MaxNumberOfCosmicClusters = cms.uint32(400000),
        MaxNumberOfPixelClusters = cms.uint32(40000),
        PixelClusterCollectionLabel = cms.InputTag("siPixelClusters"),
        cut = cms.string('strip < 400000 && pixel < 40000 && (strip < 50000 + 10*pixel) && (pixel < 5000 + 0.1*strip)'),
        doClusterCheck = cms.bool(False)
    ),
    OrderedHitsFactoryPSet = cms.PSet(
        ComponentName = cms.string('StandardHitTripletGenerator'),
        GeneratorPSet = cms.PSet(
            ComponentName = cms.string('PixelTripletHLTGenerator'),
            SeedComparitorPSet = cms.PSet(
                ComponentName = cms.string('none')
            ),
            extraHitRPhitolerance = cms.double(0.032),
            extraHitRZtolerance = cms.double(0.037),
            maxElement = cms.uint32(1000000),
            phiPreFiltering = cms.double(0.3),
            useBending = cms.bool(True),
            useFixedPreFiltering = cms.bool(False),
            useMultScattering = cms.bool(True)
        ),
        SeedingLayers = cms.InputTag("PixelLayerTriplets"),
        maxElement = cms.uint32(1000000)
    ),
    RegionFactoryPSet = cms.PSet(
        ComponentName = cms.string('GlobalRegionProducerFromBeamSpot'),
        RegionPSet = cms.PSet(
            beamSpot = cms.InputTag("offlineBeamSpot"),
            originHalfLength = cms.double(21.2),
            originRadius = cms.double(0.2),
            precise = cms.bool(True),
            ptMin = cms.double(0.9),
            useMultipleScattering = cms.bool(False)
        )
    ),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('none')
    ),
    SeedCreatorPSet = cms.PSet(
        ComponentName = cms.string('SeedFromConsecutiveHitsCreator'),
        MinOneOverPtError = cms.double(1.0),
        OriginTransverseErrorMultiplier = cms.double(1.0),
        SeedMomentumForBOFF = cms.double(5.0),
        TTRHBuilder = cms.string('WithTrackAngle'),
        forceKinematicWithRegionDirection = cms.bool(False),
        magneticField = cms.string('ParabolicMf'),
        propagator = cms.string('PropagatorWithMaterialParabolicMf')
    )
)


process.globalTrackingRegionWithVertices = cms.EDProducer("GlobalTrackingRegionWithVerticesEDProducer",
    RegionPSet = cms.PSet(
        VertexCollection = cms.InputTag("firstStepPrimaryVertices"),
        beamSpot = cms.InputTag("offlineBeamSpot"),
        fixedError = cms.double(0.2),
        halfLengthScaling4BigEvts = cms.bool(False),
        maxNVertices = cms.int32(-1),
        maxPtMin = cms.double(1000),
        minHalfLength = cms.double(0),
        minOriginR = cms.double(0),
        nSigmaZ = cms.double(4),
        originRScaling4BigEvts = cms.bool(False),
        originRadius = cms.double(0.2),
        pixelClustersForScaling = cms.InputTag("siPixelClusters"),
        precise = cms.bool(True),
        ptMin = cms.double(0.9),
        ptMinScaling4BigEvts = cms.bool(False),
        scalingEndNPix = cms.double(1),
        scalingStartNPix = cms.double(0),
        sigmaZVertex = cms.double(3),
        useFakeVertices = cms.bool(False),
        useFixedError = cms.bool(True),
        useFoundVertices = cms.bool(True),
        useMultipleScattering = cms.bool(False)
    ),
    mightGet = cms.optional.untracked.vstring
)


process.highPtTripletStep = cms.EDProducer("TrackTfClassifier",
    beamspot = cms.InputTag("offlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mightGet = cms.optional.untracked.vstring,
    mva = cms.PSet(
        batchSize = cms.int32(16),
        tfDnnLabel = cms.string('trackSelectionTf')
    ),
    qualityCuts = cms.vdouble(0.47, 0.55, 0.62),
    src = cms.InputTag("highPtTripletStepTracks"),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.highPtTripletStepClusters = cms.EDProducer("TrackClusterRemover",
    TrackQuality = cms.string('highPurity'),
    maxChi2 = cms.double(9.0),
    mightGet = cms.optional.untracked.vstring,
    minNumberOfLayersWithMeasBeforeFiltering = cms.int32(0),
    oldClusterRemovalInfo = cms.InputTag("lowPtQuadStepClusters"),
    overrideTrkQuals = cms.InputTag(""),
    pixelClusters = cms.InputTag("siPixelClusters"),
    stripClusters = cms.InputTag("siStripClusters"),
    trackClassifier = cms.InputTag("lowPtQuadStep","QualityMasks"),
    trajectories = cms.InputTag("lowPtQuadStepTracks")
)


process.highPtTripletStepHitDoublets = cms.EDProducer("HitPairEDProducer",
    clusterCheck = cms.InputTag("trackerClusterCheck"),
    layerPairs = cms.vuint32(0, 1),
    maxElement = cms.uint32(50000000),
    maxElementTotal = cms.uint32(50000000),
    mightGet = cms.optional.untracked.vstring,
    produceIntermediateHitDoublets = cms.bool(True),
    produceSeedingHitSets = cms.bool(False),
    seedingLayers = cms.InputTag("highPtTripletStepSeedLayers"),
    trackingRegions = cms.InputTag("highPtTripletStepTrackingRegions"),
    trackingRegionsSeedingLayers = cms.InputTag("")
)


process.highPtTripletStepHitTriplets = cms.EDProducer("CAHitTripletEDProducer",
    CAHardPtCut = cms.double(0.3),
    CAPhiCut = cms.double(0.07),
    CAPhiCut_byTriplets = cms.VPSet(cms.PSet(
        cut = cms.double(-1),
        seedingLayers = cms.string('')
    )),
    CAThetaCut = cms.double(0.004),
    CAThetaCut_byTriplets = cms.VPSet(cms.PSet(
        cut = cms.double(-1),
        seedingLayers = cms.string('')
    )),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('LowPtClusterShapeSeedComparitor'),
        clusterShapeCacheSrc = cms.InputTag("siPixelClusterShapeCache"),
        clusterShapeHitFilter = cms.string('ClusterShapeHitFilter')
    ),
    doublets = cms.InputTag("highPtTripletStepHitDoublets"),
    extraHitRPhitolerance = cms.double(0.032),
    maxChi2 = cms.PSet(
        enabled = cms.bool(True),
        pt1 = cms.double(0.8),
        pt2 = cms.double(8),
        value1 = cms.double(100),
        value2 = cms.double(6)
    ),
    mightGet = cms.optional.untracked.vstring,
    useBendingCorrection = cms.bool(True)
)


process.highPtTripletStepMasks = cms.EDProducer("FastTrackerRecHitMaskProducer",
    TrackQuality = cms.string('highPurity'),
    minNumberOfLayersWithMeasBeforeFiltering = cms.int32(0),
    minNumberOfLayersWithMeasBeforeFiltering_ = cms.int32(0),
    oldHitRemovalInfo = cms.InputTag("lowPtQuadStepMasks"),
    recHits = cms.InputTag("fastTrackerRecHits"),
    trackClassifier = cms.InputTag("lowPtQuadStep","QualityMasks"),
    trackQuality = cms.string('hightPurity'),
    trajectories = cms.InputTag("lowPtQuadStepTracks")
)


process.highPtTripletStepSeedClusterMask = cms.EDProducer("SeedClusterRemover",
    Common = cms.PSet(
        maxChi2 = cms.double(9.0)
    ),
    oldClusterRemovalInfo = cms.InputTag("initialStepSeedClusterMask"),
    pixelClusters = cms.InputTag("siPixelClusters"),
    stripClusters = cms.InputTag("siStripClusters"),
    trajectories = cms.InputTag("highPtTripletStepSeeds")
)


process.highPtTripletStepSeedLayers = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(
        HitProducer = cms.string('siPixelRecHits'),
        TTRHBuilder = cms.string('WithTrackAngle'),
        skipClusters = cms.InputTag("highPtTripletStepClusters")
    ),
    FPix = cms.PSet(
        HitProducer = cms.string('siPixelRecHits'),
        TTRHBuilder = cms.string('WithTrackAngle'),
        skipClusters = cms.InputTag("highPtTripletStepClusters")
    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(

    ),
    TIB = cms.PSet(

    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(

    ),
    layerList = cms.vstring(
        'BPix1+BPix2+BPix3',
        'BPix2+BPix3+BPix4',
        'BPix1+BPix3+BPix4',
        'BPix1+BPix2+BPix4',
        'BPix2+BPix3+FPix1_pos',
        'BPix2+BPix3+FPix1_neg',
        'BPix1+BPix2+FPix1_pos',
        'BPix1+BPix2+FPix1_neg',
        'BPix1+BPix3+FPix1_pos',
        'BPix1+BPix3+FPix1_neg',
        'BPix2+FPix1_pos+FPix2_pos',
        'BPix2+FPix1_neg+FPix2_neg',
        'BPix1+FPix1_pos+FPix2_pos',
        'BPix1+FPix1_neg+FPix2_neg',
        'BPix1+BPix2+FPix2_pos',
        'BPix1+BPix2+FPix2_neg',
        'FPix1_pos+FPix2_pos+FPix3_pos',
        'FPix1_neg+FPix2_neg+FPix3_neg',
        'BPix1+FPix2_pos+FPix3_pos',
        'BPix1+FPix2_neg+FPix3_neg',
        'BPix1+FPix1_pos+FPix3_pos',
        'BPix1+FPix1_neg+FPix3_neg'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.highPtTripletStepSeeds = cms.EDProducer("SeedCreatorFromRegionConsecutiveHitsEDProducer",
    MinOneOverPtError = cms.double(1),
    OriginTransverseErrorMultiplier = cms.double(1),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('none')
    ),
    SeedMomentumForBOFF = cms.double(5),
    TTRHBuilder = cms.string('WithTrackAngle'),
    forceKinematicWithRegionDirection = cms.bool(False),
    magneticField = cms.string('ParabolicMf'),
    mightGet = cms.optional.untracked.vstring,
    propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    seedingHitSets = cms.InputTag("highPtTripletStepHitTriplets")
)


process.highPtTripletStepSelector = cms.EDProducer("MultiTrackSelector",
    beamspot = cms.InputTag("offlineBeamSpot"),
    src = cms.InputTag("highPtTripletStepTracks"),
    trackSelectors = cms.VPSet(
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(True),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(2.0),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(0.7, 4.0),
            d0_par2 = cms.vdouble(0.6, 4.0),
            dz_par1 = cms.vdouble(0.8, 4.0),
            dz_par2 = cms.vdouble(0.6, 4.0),
            keepAllTracks = cms.bool(False),
            maxNumberLostLayers = cms.uint32(3),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(3),
            minNumberLayers = cms.uint32(3),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(0),
            nSigmaZ = cms.double(4.0),
            name = cms.string('highPtTripletStepLoose'),
            preFilterName = cms.string(''),
            qualityBit = cms.string('loose'),
            res_par = cms.vdouble(0.003, 0.002),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        ),
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(True),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(1.0),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(0.6, 4.0),
            d0_par2 = cms.vdouble(0.5, 4.0),
            dz_par1 = cms.vdouble(0.7, 4.0),
            dz_par2 = cms.vdouble(0.6, 4.0),
            keepAllTracks = cms.bool(True),
            maxNumberLostLayers = cms.uint32(2),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(3),
            minNumberLayers = cms.uint32(3),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(0),
            nSigmaZ = cms.double(4.0),
            name = cms.string('highPtTripletStepTight'),
            preFilterName = cms.string('highPtTripletStepLoose'),
            qualityBit = cms.string('tight'),
            res_par = cms.vdouble(0.003, 0.002),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        ),
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(True),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(0.8),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(0.6, 4.0),
            d0_par2 = cms.vdouble(0.45, 4.0),
            dz_par1 = cms.vdouble(0.7, 4.0),
            dz_par2 = cms.vdouble(0.55, 4.0),
            keepAllTracks = cms.bool(True),
            maxNumberLostLayers = cms.uint32(2),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(4.1),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(4),
            minNumberLayers = cms.uint32(4),
            min_eta = cms.double(-4.1),
            min_nhits = cms.uint32(4),
            nSigmaZ = cms.double(4.0),
            name = cms.string('highPtTripletStep'),
            preFilterName = cms.string('highPtTripletStepTight'),
            qualityBit = cms.string('highPurity'),
            res_par = cms.vdouble(0.003, 0.001),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        )
    ),
    useVertices = cms.bool(True),
    useVtxError = cms.bool(False),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.highPtTripletStepTrackCandidates = cms.EDProducer("MkFitOutputConverter",
    batchSize = cms.int32(16),
    candMVASel = cms.bool(True),
    candWP = cms.double(-0.3),
    doErrorRescale = cms.bool(True),
    mightGet = cms.optional.untracked.vstring,
    mkFitEventOfHits = cms.InputTag("mkFitEventOfHits"),
    mkFitPixelHits = cms.InputTag("mkFitSiPixelHits"),
    mkFitSeeds = cms.InputTag("highPtTripletStepTrackCandidatesMkFitSeeds"),
    mkFitStripHits = cms.InputTag("mkFitSiStripHits"),
    propagatorAlong = cms.ESInputTag("","PropagatorWithMaterial"),
    propagatorOpposite = cms.ESInputTag("","PropagatorWithMaterialOpposite"),
    qualityMaxInvPt = cms.double(100),
    qualityMaxPosErr = cms.double(100),
    qualityMaxR = cms.double(120),
    qualityMaxZ = cms.double(280),
    qualityMinTheta = cms.double(0.01),
    qualitySignPt = cms.bool(True),
    seeds = cms.InputTag("highPtTripletStepSeeds"),
    tfDnnLabel = cms.string('trackSelectionTf'),
    tracks = cms.InputTag("highPtTripletStepTrackCandidatesMkFit"),
    ttrhBuilder = cms.ESInputTag("","WithTrackAngle")
)


process.highPtTripletStepTrackCandidatesMkFit = cms.EDProducer("MkFitProducer",
    backwardFitInCMSSW = cms.bool(False),
    buildingRoutine = cms.string('cloneEngine'),
    clustersToSkip = cms.InputTag("highPtTripletStepClusters"),
    config = cms.ESInputTag("","highPtTripletStepTrackCandidatesMkFitConfig"),
    eventOfHits = cms.InputTag("mkFitEventOfHits"),
    limitConcurrency = cms.untracked.bool(False),
    mightGet = cms.optional.untracked.vstring,
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutLoose')
    ),
    mkFitSilent = cms.untracked.bool(True),
    pixelHits = cms.InputTag("mkFitSiPixelHits"),
    removeDuplicates = cms.bool(True),
    seedCleaning = cms.bool(True),
    seeds = cms.InputTag("highPtTripletStepTrackCandidatesMkFitSeeds"),
    stripHits = cms.InputTag("mkFitSiStripHits")
)


process.highPtTripletStepTrackCandidatesMkFitSeeds = cms.EDProducer("MkFitSeedConverter",
    maxNSeeds = cms.uint32(500000),
    mightGet = cms.optional.untracked.vstring,
    seeds = cms.InputTag("highPtTripletStepSeeds"),
    ttrhBuilder = cms.ESInputTag("","WithTrackAngle")
)


process.highPtTripletStepTrackingRegions = cms.EDProducer("GlobalTrackingRegionFromBeamSpotEDProducer",
    RegionPSet = cms.PSet(
        beamSpot = cms.InputTag("offlineBeamSpot"),
        nSigmaZ = cms.double(4),
        originHalfLength = cms.double(0),
        originRadius = cms.double(0.02),
        precise = cms.bool(True),
        ptMin = cms.double(0.55),
        useMultipleScattering = cms.bool(False)
    ),
    mightGet = cms.optional.untracked.vstring
)


process.highPtTripletStepTracks = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('highPtTripletStep'),
    Fitter = cms.string('FlexibleKFFittingSmoother'),
    GeometricInnerState = cms.bool(False),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    Propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    TTRHBuilder = cms.string('WithAngleAndTemplateWithoutProbQ'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("highPtTripletStepTrackCandidates"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(True)
)


process.hpsPFTauBasicDiscriminators = cms.EDProducer("PFRecoTauDiscriminationByIsolationContainer",
    IDWPdefinitions = cms.VPSet(
        cms.PSet(
            IDname = cms.string('ByLooseCombinedIsolationDBSumPtCorr3Hits'),
            maximumAbsoluteValues = cms.vdouble(2.5, 1000000000.0),
            maximumRelativeValues = cms.vdouble(-1.0, 0.1),
            referenceRawIDNames = cms.vstring(
                'ByRawCombinedIsolationDBSumPtCorr3Hits',
                'PhotonPtSumOutsideSignalCone'
            )
        ),
        cms.PSet(
            IDname = cms.string('ByMediumCombinedIsolationDBSumPtCorr3Hits'),
            maximumAbsoluteValues = cms.vdouble(1.5, 1000000000.0),
            maximumRelativeValues = cms.vdouble(-1.0, 0.1),
            referenceRawIDNames = cms.vstring(
                'ByRawCombinedIsolationDBSumPtCorr3Hits',
                'PhotonPtSumOutsideSignalCone'
            )
        ),
        cms.PSet(
            IDname = cms.string('ByTightCombinedIsolationDBSumPtCorr3Hits'),
            maximumAbsoluteValues = cms.vdouble(0.8, 1000000000.0),
            maximumRelativeValues = cms.vdouble(-1.0, 0.1),
            referenceRawIDNames = cms.vstring(
                'ByRawCombinedIsolationDBSumPtCorr3Hits',
                'PhotonPtSumOutsideSignalCone'
            )
        ),
        cms.PSet(
            IDname = cms.string('ByLooseChargedIsolation'),
            maximumAbsoluteValues = cms.vdouble(2.5),
            referenceRawIDNames = cms.vstring('ChargedIsoPtSum')
        ),
        cms.PSet(
            IDname = cms.string('ByPhotonPtSumOutsideSignalCone'),
            maximumRelativeValues = cms.vdouble(0.1),
            referenceRawIDNames = cms.vstring('PhotonPtSumOutsideSignalCone')
        )
    ),
    IDdefinitions = cms.VPSet(
        cms.PSet(
            ApplyDiscriminationByTrackerIsolation = cms.bool(True),
            IDname = cms.string('ChargedIsoPtSum'),
            storeRawSumPt = cms.bool(True)
        ),
        cms.PSet(
            ApplyDiscriminationByECALIsolation = cms.bool(True),
            IDname = cms.string('NeutralIsoPtSum'),
            storeRawSumPt = cms.bool(True)
        ),
        cms.PSet(
            ApplyDiscriminationByWeightedECALIsolation = cms.bool(True),
            IDname = cms.string('NeutralIsoPtSumWeight'),
            UseAllPFCandsForWeights = cms.bool(True),
            storeRawSumPt = cms.bool(True)
        ),
        cms.PSet(
            IDname = cms.string('TauFootprintCorrection'),
            storeRawFootprintCorrection = cms.bool(True)
        ),
        cms.PSet(
            IDname = cms.string('PhotonPtSumOutsideSignalCone'),
            storeRawPhotonSumPt_outsideSignalCone = cms.bool(True)
        ),
        cms.PSet(
            IDname = cms.string('PUcorrPtSum'),
            applyDeltaBetaCorrection = cms.bool(True),
            storeRawPUsumPt = cms.bool(True)
        ),
        cms.PSet(
            ApplyDiscriminationByECALIsolation = cms.bool(True),
            ApplyDiscriminationByTrackerIsolation = cms.bool(True),
            IDname = cms.string('ByRawCombinedIsolationDBSumPtCorr3Hits'),
            applyDeltaBetaCorrection = cms.bool(True),
            storeRawSumPt = cms.bool(True)
        )
    ),
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        decayMode = cms.PSet(
            Producer = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs"),
            cut = cms.double(0.5)
        )
    ),
    WeightECALIsolation = cms.double(1),
    applyFootprintCorrection = cms.bool(False),
    applyRhoCorrection = cms.bool(False),
    customOuterCone = cms.double(0.5),
    deltaBetaFactor = cms.string('0.2000'),
    deltaBetaPUTrackPtCutOverride = cms.bool(True),
    deltaBetaPUTrackPtCutOverride_val = cms.double(0.5),
    footprintCorrections = cms.VPSet(
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 0')
        ),
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 1 || decayMode() = 2')
        ),
        cms.PSet(
            offset = cms.string('2.7'),
            selection = cms.string('decayMode() = 5')
        ),
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 6')
        ),
        cms.PSet(
            offset = cms.string('max(2.0, 0.22*pt() - 2.0)'),
            selection = cms.string('decayMode() = 10')
        )
    ),
    isoConeSizeForDeltaBeta = cms.double(0.8),
    mightGet = cms.optional.untracked.vstring,
    minTauPtForNoIso = cms.double(-99),
    particleFlowSrc = cms.InputTag("particleFlow"),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.2),
            maxDeltaZToLeadTrack = cms.double(-1.0),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.03),
            minGammaEt = cms.double(1.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        leadingTrkOrPFCandOption = cms.string('leadPFCand'),
        primaryVertexSrc = cms.InputTag("offlinePrimaryVertices"),
        pvFindingAlgo = cms.string('closestInDeltaZ'),
        recoverLeadingTrk = cms.bool(False),
        signalQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.4),
            maxDeltaZToLeadTrack = cms.double(-1.0),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(1.0),
            minNeutralHadronEt = cms.double(30.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        vertexTrackFiltering = cms.bool(False),
        vxAssocQualityCuts = cms.PSet(
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(1.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        )
    ),
    rhoConeSize = cms.double(0.5),
    rhoProducer = cms.InputTag("fixedGridRhoFastjetAll"),
    rhoUEOffsetCorrection = cms.double(1),
    verbosity = cms.int32(0),
    vertexSrc = cms.InputTag("offlinePrimaryVertices")
)


process.hpsPFTauBasicDiscriminatorsdR03 = cms.EDProducer("PFRecoTauDiscriminationByIsolationContainer",
    IDWPdefinitions = cms.VPSet(
        cms.PSet(
            IDname = cms.string('ByLooseCombinedIsolationDBSumPtCorr3HitsdR03'),
            maximumAbsoluteValues = cms.vdouble(2.5, 1000000000.0),
            maximumRelativeValues = cms.vdouble(-1.0, 0.1),
            referenceRawIDNames = cms.vstring(
                'ByRawCombinedIsolationDBSumPtCorr3HitsdR03',
                'PhotonPtSumOutsideSignalConedR03'
            )
        ),
        cms.PSet(
            IDname = cms.string('ByMediumCombinedIsolationDBSumPtCorr3HitsdR03'),
            maximumAbsoluteValues = cms.vdouble(1.5, 1000000000.0),
            maximumRelativeValues = cms.vdouble(-1.0, 0.1),
            referenceRawIDNames = cms.vstring(
                'ByRawCombinedIsolationDBSumPtCorr3HitsdR03',
                'PhotonPtSumOutsideSignalConedR03'
            )
        ),
        cms.PSet(
            IDname = cms.string('ByTightCombinedIsolationDBSumPtCorr3HitsdR03'),
            maximumAbsoluteValues = cms.vdouble(0.8, 1000000000.0),
            maximumRelativeValues = cms.vdouble(-1.0, 0.1),
            referenceRawIDNames = cms.vstring(
                'ByRawCombinedIsolationDBSumPtCorr3HitsdR03',
                'PhotonPtSumOutsideSignalConedR03'
            )
        )
    ),
    IDdefinitions = cms.VPSet(
        cms.PSet(
            ApplyDiscriminationByTrackerIsolation = cms.bool(True),
            IDname = cms.string('ChargedIsoPtSumdR03'),
            storeRawSumPt = cms.bool(True)
        ),
        cms.PSet(
            ApplyDiscriminationByECALIsolation = cms.bool(True),
            IDname = cms.string('NeutralIsoPtSumdR03'),
            storeRawSumPt = cms.bool(True)
        ),
        cms.PSet(
            ApplyDiscriminationByWeightedECALIsolation = cms.bool(True),
            IDname = cms.string('NeutralIsoPtSumWeightdR03'),
            UseAllPFCandsForWeights = cms.bool(True),
            storeRawSumPt = cms.bool(True)
        ),
        cms.PSet(
            IDname = cms.string('TauFootprintCorrectiondR03'),
            storeRawFootprintCorrection = cms.bool(True)
        ),
        cms.PSet(
            IDname = cms.string('PhotonPtSumOutsideSignalConedR03'),
            storeRawPhotonSumPt_outsideSignalCone = cms.bool(True)
        ),
        cms.PSet(
            IDname = cms.string('PUcorrPtSumdR03'),
            applyDeltaBetaCorrection = cms.bool(True),
            storeRawPUsumPt = cms.bool(True)
        ),
        cms.PSet(
            ApplyDiscriminationByECALIsolation = cms.bool(True),
            ApplyDiscriminationByTrackerIsolation = cms.bool(True),
            IDname = cms.string('ByRawCombinedIsolationDBSumPtCorr3HitsdR03'),
            applyDeltaBetaCorrection = cms.bool(True),
            storeRawSumPt = cms.bool(True)
        )
    ),
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        decayMode = cms.PSet(
            Producer = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs"),
            cut = cms.double(0.5)
        )
    ),
    WeightECALIsolation = cms.double(1),
    applyFootprintCorrection = cms.bool(False),
    applyRhoCorrection = cms.bool(False),
    customOuterCone = cms.double(0.3),
    deltaBetaFactor = cms.string('0.0720'),
    deltaBetaPUTrackPtCutOverride = cms.bool(True),
    deltaBetaPUTrackPtCutOverride_val = cms.double(0.5),
    footprintCorrections = cms.VPSet(
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 0')
        ),
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 1 || decayMode() = 2')
        ),
        cms.PSet(
            offset = cms.string('2.7'),
            selection = cms.string('decayMode() = 5')
        ),
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 6')
        ),
        cms.PSet(
            offset = cms.string('max(2.0, 0.22*pt() - 2.0)'),
            selection = cms.string('decayMode() = 10')
        )
    ),
    isoConeSizeForDeltaBeta = cms.double(0.8),
    mightGet = cms.optional.untracked.vstring,
    minTauPtForNoIso = cms.double(-99),
    particleFlowSrc = cms.InputTag("particleFlow"),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.2),
            maxDeltaZToLeadTrack = cms.double(-1.0),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.03),
            minGammaEt = cms.double(1.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        leadingTrkOrPFCandOption = cms.string('leadPFCand'),
        primaryVertexSrc = cms.InputTag("offlinePrimaryVertices"),
        pvFindingAlgo = cms.string('closestInDeltaZ'),
        recoverLeadingTrk = cms.bool(False),
        signalQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.4),
            maxDeltaZToLeadTrack = cms.double(-1.0),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(1.0),
            minNeutralHadronEt = cms.double(30.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        vertexTrackFiltering = cms.bool(False),
        vxAssocQualityCuts = cms.PSet(
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(1.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        )
    ),
    rhoConeSize = cms.double(0.5),
    rhoProducer = cms.InputTag("fixedGridRhoFastjetAll"),
    rhoUEOffsetCorrection = cms.double(1),
    verbosity = cms.int32(0),
    vertexSrc = cms.InputTag("offlinePrimaryVertices")
)


process.hpsPFTauDiscriminationByDeadECALElectronRejection = cms.EDProducer("PFRecoTauDiscriminationAgainstElectronDeadECAL",
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        decayMode = cms.PSet(
            Producer = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs"),
            cut = cms.double(0.5)
        )
    ),
    dR = cms.double(0.08),
    extrapolateToECalEntrance = cms.bool(True),
    mightGet = cms.optional.untracked.vstring,
    minStatus = cms.uint32(12),
    verbosity = cms.int32(0)
)


process.hpsPFTauDiscriminationByDecayModeFinding = cms.EDProducer("PFRecoTauDiscriminationByHPSSelection",
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and')
    ),
    decayModes = cms.VPSet(
        cms.PSet(
            applyBendCorrection = cms.PSet(
                eta = cms.bool(True),
                mass = cms.bool(True),
                phi = cms.bool(True)
            ),
            assumeStripMass = cms.double(-1.0),
            maxMass = cms.string('1.'),
            maxPi0Mass = cms.double(1000000000.0),
            minMass = cms.double(-1000.0),
            minPi0Mass = cms.double(-1000.0),
            nCharged = cms.uint32(1),
            nChargedPFCandsMin = cms.uint32(1),
            nPiZeros = cms.uint32(0),
            nTracksMin = cms.uint32(1)
        ),
        cms.PSet(
            applyBendCorrection = cms.PSet(
                eta = cms.bool(True),
                mass = cms.bool(True),
                phi = cms.bool(True)
            ),
            assumeStripMass = cms.double(0.1349),
            maxMass = cms.string('max(1.3, min(1.3*sqrt(pt/100.), 4.2))'),
            maxPi0Mass = cms.double(1000000000.0),
            minMass = cms.double(0.3),
            minPi0Mass = cms.double(-1000.0),
            nCharged = cms.uint32(1),
            nChargedPFCandsMin = cms.uint32(1),
            nPiZeros = cms.uint32(1),
            nTracksMin = cms.uint32(1)
        ),
        cms.PSet(
            applyBendCorrection = cms.PSet(
                eta = cms.bool(True),
                mass = cms.bool(True),
                phi = cms.bool(True)
            ),
            assumeStripMass = cms.double(0.0),
            maxMass = cms.string('max(1.2, min(1.2*sqrt(pt/100.), 4.0))'),
            maxPi0Mass = cms.double(0.2),
            minMass = cms.double(0.4),
            minPi0Mass = cms.double(0.05),
            nCharged = cms.uint32(1),
            nChargedPFCandsMin = cms.uint32(1),
            nPiZeros = cms.uint32(2),
            nTracksMin = cms.uint32(1)
        ),
        cms.PSet(
            applyBendCorrection = cms.PSet(
                eta = cms.bool(False),
                mass = cms.bool(False),
                phi = cms.bool(False)
            ),
            assumeStripMass = cms.double(-1.0),
            maxMass = cms.string('1.5'),
            maxPi0Mass = cms.double(1000000000.0),
            minMass = cms.double(0.8),
            minPi0Mass = cms.double(-1000.0),
            nCharged = cms.uint32(3),
            nChargedPFCandsMin = cms.uint32(1),
            nPiZeros = cms.uint32(0),
            nTracksMin = cms.uint32(2)
        )
    ),
    matchingCone = cms.double(0.5),
    minPixelHits = cms.int32(1),
    minTauPt = cms.double(0.0),
    requireTauChargedHadronsToBeChargedPFCands = cms.bool(True),
    verbosity = cms.int32(0)
)


process.hpsPFTauDiscriminationByDecayModeFindingNewDMs = cms.EDProducer("PFRecoTauDiscriminationByHPSSelection",
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and')
    ),
    decayModes = cms.VPSet(
        cms.PSet(
            applyBendCorrection = cms.PSet(
                eta = cms.bool(True),
                mass = cms.bool(True),
                phi = cms.bool(True)
            ),
            assumeStripMass = cms.double(-1.0),
            maxMass = cms.string('1.'),
            maxPi0Mass = cms.double(1000000000.0),
            minMass = cms.double(-1000.0),
            minPi0Mass = cms.double(-1000.0),
            nCharged = cms.uint32(1),
            nChargedPFCandsMin = cms.uint32(1),
            nPiZeros = cms.uint32(0),
            nTracksMin = cms.uint32(1)
        ),
        cms.PSet(
            applyBendCorrection = cms.PSet(
                eta = cms.bool(True),
                mass = cms.bool(True),
                phi = cms.bool(True)
            ),
            assumeStripMass = cms.double(0.1349),
            maxMass = cms.string('max(1.3, min(1.3*sqrt(pt/100.), 4.2))'),
            maxPi0Mass = cms.double(1000000000.0),
            minMass = cms.double(0.3),
            minPi0Mass = cms.double(-1000.0),
            nCharged = cms.uint32(1),
            nChargedPFCandsMin = cms.uint32(1),
            nPiZeros = cms.uint32(1),
            nTracksMin = cms.uint32(1)
        ),
        cms.PSet(
            applyBendCorrection = cms.PSet(
                eta = cms.bool(True),
                mass = cms.bool(True),
                phi = cms.bool(True)
            ),
            assumeStripMass = cms.double(0.0),
            maxMass = cms.string('max(1.2, min(1.2*sqrt(pt/100.), 4.0))'),
            maxPi0Mass = cms.double(0.2),
            minMass = cms.double(0.4),
            minPi0Mass = cms.double(0.05),
            nCharged = cms.uint32(1),
            nChargedPFCandsMin = cms.uint32(1),
            nPiZeros = cms.uint32(2),
            nTracksMin = cms.uint32(1)
        ),
        cms.PSet(
            applyBendCorrection = cms.PSet(
                eta = cms.bool(False),
                mass = cms.bool(False),
                phi = cms.bool(False)
            ),
            assumeStripMass = cms.double(-1.0),
            maxMass = cms.string('1.2'),
            maxPi0Mass = cms.double(1000000000.0),
            minMass = cms.double(0.0),
            minPi0Mass = cms.double(-1000.0),
            nCharged = cms.uint32(2),
            nChargedPFCandsMin = cms.uint32(1),
            nPiZeros = cms.uint32(0),
            nTracksMin = cms.uint32(2)
        ),
        cms.PSet(
            applyBendCorrection = cms.PSet(
                eta = cms.bool(False),
                mass = cms.bool(False),
                phi = cms.bool(False)
            ),
            assumeStripMass = cms.double(-1.0),
            maxMass = cms.string('max(1.2, min(1.2*sqrt(pt/100.), 4.0))'),
            maxPi0Mass = cms.double(1000000000.0),
            minMass = cms.double(0.0),
            minPi0Mass = cms.double(-1000.0),
            nCharged = cms.uint32(2),
            nChargedPFCandsMin = cms.uint32(1),
            nPiZeros = cms.uint32(1),
            nTracksMin = cms.uint32(2)
        ),
        cms.PSet(
            applyBendCorrection = cms.PSet(
                eta = cms.bool(False),
                mass = cms.bool(False),
                phi = cms.bool(False)
            ),
            assumeStripMass = cms.double(-1.0),
            maxMass = cms.string('1.5'),
            maxPi0Mass = cms.double(1000000000.0),
            minMass = cms.double(0.8),
            minPi0Mass = cms.double(-1000.0),
            nCharged = cms.uint32(3),
            nChargedPFCandsMin = cms.uint32(1),
            nPiZeros = cms.uint32(0),
            nTracksMin = cms.uint32(2)
        ),
        cms.PSet(
            applyBendCorrection = cms.PSet(
                eta = cms.bool(False),
                mass = cms.bool(False),
                phi = cms.bool(False)
            ),
            assumeStripMass = cms.double(-1.0),
            maxMass = cms.string('1.6'),
            maxPi0Mass = cms.double(1000000000.0),
            minMass = cms.double(0.9),
            minPi0Mass = cms.double(-1000.0),
            nCharged = cms.uint32(3),
            nChargedPFCandsMin = cms.uint32(1),
            nPiZeros = cms.uint32(1),
            nTracksMin = cms.uint32(2)
        )
    ),
    matchingCone = cms.double(0.5),
    minPixelHits = cms.int32(1),
    minTauPt = cms.double(0.0),
    requireTauChargedHadronsToBeChargedPFCands = cms.bool(False),
    verbosity = cms.int32(0)
)


process.hpsPFTauDiscriminationByDecayModeFindingOldDMs = cms.EDProducer("PFRecoTauDiscriminationByHPSSelection",
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and')
    ),
    decayModes = cms.VPSet(
        cms.PSet(
            applyBendCorrection = cms.PSet(
                eta = cms.bool(True),
                mass = cms.bool(True),
                phi = cms.bool(True)
            ),
            assumeStripMass = cms.double(-1.0),
            maxMass = cms.string('1.'),
            maxPi0Mass = cms.double(1000000000.0),
            minMass = cms.double(-1000.0),
            minPi0Mass = cms.double(-1000.0),
            nCharged = cms.uint32(1),
            nChargedPFCandsMin = cms.uint32(1),
            nPiZeros = cms.uint32(0),
            nTracksMin = cms.uint32(1)
        ),
        cms.PSet(
            applyBendCorrection = cms.PSet(
                eta = cms.bool(True),
                mass = cms.bool(True),
                phi = cms.bool(True)
            ),
            assumeStripMass = cms.double(0.1349),
            maxMass = cms.string('max(1.3, min(1.3*sqrt(pt/100.), 4.2))'),
            maxPi0Mass = cms.double(1000000000.0),
            minMass = cms.double(0.3),
            minPi0Mass = cms.double(-1000.0),
            nCharged = cms.uint32(1),
            nChargedPFCandsMin = cms.uint32(1),
            nPiZeros = cms.uint32(1),
            nTracksMin = cms.uint32(1)
        ),
        cms.PSet(
            applyBendCorrection = cms.PSet(
                eta = cms.bool(True),
                mass = cms.bool(True),
                phi = cms.bool(True)
            ),
            assumeStripMass = cms.double(0.0),
            maxMass = cms.string('max(1.2, min(1.2*sqrt(pt/100.), 4.0))'),
            maxPi0Mass = cms.double(0.2),
            minMass = cms.double(0.4),
            minPi0Mass = cms.double(0.05),
            nCharged = cms.uint32(1),
            nChargedPFCandsMin = cms.uint32(1),
            nPiZeros = cms.uint32(2),
            nTracksMin = cms.uint32(1)
        ),
        cms.PSet(
            applyBendCorrection = cms.PSet(
                eta = cms.bool(False),
                mass = cms.bool(False),
                phi = cms.bool(False)
            ),
            assumeStripMass = cms.double(-1.0),
            maxMass = cms.string('1.5'),
            maxPi0Mass = cms.double(1000000000.0),
            minMass = cms.double(0.8),
            minPi0Mass = cms.double(-1000.0),
            nCharged = cms.uint32(3),
            nChargedPFCandsMin = cms.uint32(1),
            nPiZeros = cms.uint32(0),
            nTracksMin = cms.uint32(2)
        )
    ),
    matchingCone = cms.double(0.5),
    minPixelHits = cms.int32(1),
    minTauPt = cms.double(0.0),
    requireTauChargedHadronsToBeChargedPFCands = cms.bool(True),
    verbosity = cms.int32(0)
)


process.hpsPFTauDiscriminationByIsolationMVArun2v1DBdR03oldDMwLT = cms.EDProducer("RecoTauDiscriminantCutMultiplexer",
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        decayMode = cms.PSet(
            Producer = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs"),
            cut = cms.double(0.5)
        )
    ),
    inputFileName = cms.FileInPath('RecoTauTag/RecoTau/data/emptyMVAinputFile'),
    loadMVAfromDB = cms.bool(True),
    mapping = cms.VPSet(cms.PSet(
        category = cms.uint32(0),
        cut = cms.string('RecoTauTag_tauIdMVAIsoDBoldDMdR0p3wLT'),
        variable = cms.string('pt')
    )),
    mightGet = cms.optional.untracked.vstring,
    mvaOutput_normalization = cms.string('RecoTauTag_tauIdMVAIsoDBoldDMdR0p3wLT_mvaOutput_normalization'),
    rawValues = cms.vstring('discriminator'),
    toMultiplex = cms.InputTag("hpsPFTauDiscriminationByIsolationMVArun2v1DBdR03oldDMwLTraw"),
    verbosity = cms.int32(0),
    workingPoints = cms.vstring(
        '_VVLoose',
        '_VLoose',
        '_Loose',
        '_Medium',
        '_Tight',
        '_VTight',
        '_VVTight'
    )
)


process.hpsPFTauDiscriminationByIsolationMVArun2v1DBdR03oldDMwLTraw = cms.EDProducer("PFRecoTauDiscriminationByMVAIsolationRun2",
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        decayMode = cms.PSet(
            Producer = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs"),
            cut = cms.double(0.5)
        )
    ),
    inputFileName = cms.optional.FileInPath,
    inputIDNameSuffix = cms.string('dR03'),
    loadMVAfromDB = cms.bool(True),
    mightGet = cms.optional.untracked.vstring,
    mvaName = cms.string('RecoTauTag_tauIdMVAIsoDBoldDMdR0p3wLT'),
    mvaOpt = cms.string('DBoldDMwLTwGJ'),
    srcBasicTauDiscriminators = cms.InputTag("hpsPFTauBasicDiscriminatorsdR03"),
    srcTauTransverseImpactParameters = cms.InputTag("hpsPFTauTransverseImpactParameters"),
    verbosity = cms.int32(0)
)


process.hpsPFTauDiscriminationByIsolationMVArun2v1DBnewDMwLT = cms.EDProducer("RecoTauDiscriminantCutMultiplexer",
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        decayMode = cms.PSet(
            Producer = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs"),
            cut = cms.double(0.5)
        )
    ),
    inputFileName = cms.FileInPath('RecoTauTag/RecoTau/data/emptyMVAinputFile'),
    loadMVAfromDB = cms.bool(True),
    mapping = cms.VPSet(cms.PSet(
        category = cms.uint32(0),
        cut = cms.string('RecoTauTag_tauIdMVAIsoDBnewDMwLT'),
        variable = cms.string('pt')
    )),
    mightGet = cms.optional.untracked.vstring,
    mvaOutput_normalization = cms.string('RecoTauTag_tauIdMVAIsoDBnewDMwLT_mvaOutput_normalization'),
    rawValues = cms.vstring('discriminator'),
    toMultiplex = cms.InputTag("hpsPFTauDiscriminationByIsolationMVArun2v1DBnewDMwLTraw"),
    verbosity = cms.int32(0),
    workingPoints = cms.vstring(
        '_VVLoose',
        '_VLoose',
        '_Loose',
        '_Medium',
        '_Tight',
        '_VTight',
        '_VVTight'
    )
)


process.hpsPFTauDiscriminationByIsolationMVArun2v1DBnewDMwLTraw = cms.EDProducer("PFRecoTauDiscriminationByMVAIsolationRun2",
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        decayMode = cms.PSet(
            Producer = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs"),
            cut = cms.double(0.5)
        )
    ),
    inputFileName = cms.optional.FileInPath,
    inputIDNameSuffix = cms.string(''),
    loadMVAfromDB = cms.bool(True),
    mightGet = cms.optional.untracked.vstring,
    mvaName = cms.string('RecoTauTag_tauIdMVAIsoDBnewDMwLT'),
    mvaOpt = cms.string('DBnewDMwLTwGJ'),
    srcBasicTauDiscriminators = cms.InputTag("hpsPFTauBasicDiscriminators"),
    srcTauTransverseImpactParameters = cms.InputTag("hpsPFTauTransverseImpactParameters"),
    verbosity = cms.int32(0)
)


process.hpsPFTauDiscriminationByIsolationMVArun2v1DBoldDMwLT = cms.EDProducer("RecoTauDiscriminantCutMultiplexer",
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        decayMode = cms.PSet(
            Producer = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs"),
            cut = cms.double(0.5)
        )
    ),
    inputFileName = cms.FileInPath('RecoTauTag/RecoTau/data/emptyMVAinputFile'),
    loadMVAfromDB = cms.bool(True),
    mapping = cms.VPSet(cms.PSet(
        category = cms.uint32(0),
        cut = cms.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT'),
        variable = cms.string('pt')
    )),
    mightGet = cms.optional.untracked.vstring,
    mvaOutput_normalization = cms.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT_mvaOutput_normalization'),
    rawValues = cms.vstring('discriminator'),
    toMultiplex = cms.InputTag("hpsPFTauDiscriminationByIsolationMVArun2v1DBoldDMwLTraw"),
    verbosity = cms.int32(0),
    workingPoints = cms.vstring(
        '_VVLoose',
        '_VLoose',
        '_Loose',
        '_Medium',
        '_Tight',
        '_VTight',
        '_VVTight'
    )
)


process.hpsPFTauDiscriminationByIsolationMVArun2v1DBoldDMwLTraw = cms.EDProducer("PFRecoTauDiscriminationByMVAIsolationRun2",
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        decayMode = cms.PSet(
            Producer = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs"),
            cut = cms.double(0.5)
        )
    ),
    inputFileName = cms.optional.FileInPath,
    inputIDNameSuffix = cms.string(''),
    loadMVAfromDB = cms.bool(True),
    mightGet = cms.optional.untracked.vstring,
    mvaName = cms.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT'),
    mvaOpt = cms.string('DBoldDMwLTwGJ'),
    srcBasicTauDiscriminators = cms.InputTag("hpsPFTauBasicDiscriminators"),
    srcTauTransverseImpactParameters = cms.InputTag("hpsPFTauTransverseImpactParameters"),
    verbosity = cms.int32(0)
)


process.hpsPFTauDiscriminationByLooseElectronRejection = cms.EDProducer("PFRecoTauDiscriminationAgainstElectron",
    ApplyCut_BremCombined = cms.bool(False),
    ApplyCut_BremsRecoveryEOverPLead = cms.bool(False),
    ApplyCut_EOverPLead = cms.bool(False),
    ApplyCut_EcalCrackCut = cms.bool(False),
    ApplyCut_ElectronPreID = cms.bool(False),
    ApplyCut_ElectronPreID_2D = cms.bool(False),
    ApplyCut_EmFraction = cms.bool(False),
    ApplyCut_Hcal3x3OverPLead = cms.bool(False),
    ApplyCut_HcalMaxOverPLead = cms.bool(False),
    ApplyCut_HcalTotOverPLead = cms.bool(False),
    ApplyCut_PFElectronMVA = cms.bool(True),
    BremCombined_Fraction = cms.double(0.99),
    BremCombined_HOP = cms.double(0.1),
    BremCombined_Mass = cms.double(0.55),
    BremCombined_StripSize = cms.double(0.03),
    BremsRecoveryEOverPLead_maxValue = cms.double(1.8),
    BremsRecoveryEOverPLead_minValue = cms.double(0.8),
    EOverPLead_maxValue = cms.double(1.8),
    EOverPLead_minValue = cms.double(0.8),
    ElecPreID0_EOverPLead_maxValue = cms.double(0.95),
    ElecPreID0_HOverPLead_minValue = cms.double(0.05),
    ElecPreID1_EOverPLead_maxValue = cms.double(0.8),
    ElecPreID1_HOverPLead_minValue = cms.double(0.15),
    EmFraction_maxValue = cms.double(0.9),
    Hcal3x3OverPLead_minValue = cms.double(0.1),
    HcalMaxOverPLead_minValue = cms.double(0.1),
    HcalTotOverPLead_minValue = cms.double(0.1),
    PFElectronMVA_maxValue = cms.double(0.6),
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and')
    )
)


process.hpsPFTauDiscriminationByMVA6ElectronRejection = cms.EDProducer("RecoTauDiscriminantCutMultiplexer",
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        decayMode = cms.PSet(
            Producer = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs"),
            cut = cms.double(0.5)
        )
    ),
    inputFileName = cms.FileInPath('RecoTauTag/RecoTau/data/emptyMVAinputFile'),
    loadMVAfromDB = cms.bool(True),
    mapping = cms.VPSet(
        cms.PSet(
            category = cms.uint32(0),
            cut = cms.string('RecoTauTag_antiElectronMVA_NoEleMatch_woGwoGSF_BL'),
            variable = cms.string('pt')
        ),
        cms.PSet(
            category = cms.uint32(2),
            cut = cms.string('RecoTauTag_antiElectronMVA_NoEleMatch_wGwoGSF_BL'),
            variable = cms.string('pt')
        ),
        cms.PSet(
            category = cms.uint32(5),
            cut = cms.string('RecoTauTag_antiElectronMVA_woGwGSF_BL'),
            variable = cms.string('pt')
        ),
        cms.PSet(
            category = cms.uint32(7),
            cut = cms.string('RecoTauTag_antiElectronMVA_wGwGSF_BL'),
            variable = cms.string('pt')
        ),
        cms.PSet(
            category = cms.uint32(8),
            cut = cms.string('RecoTauTag_antiElectronMVA_NoEleMatch_woGwoGSF_EC'),
            variable = cms.string('pt')
        ),
        cms.PSet(
            category = cms.uint32(10),
            cut = cms.string('RecoTauTag_antiElectronMVA_NoEleMatch_wGwoGSF_EC'),
            variable = cms.string('pt')
        ),
        cms.PSet(
            category = cms.uint32(13),
            cut = cms.string('RecoTauTag_antiElectronMVA_woGwGSF_EC'),
            variable = cms.string('pt')
        ),
        cms.PSet(
            category = cms.uint32(15),
            cut = cms.string('RecoTauTag_antiElectronMVA_wGwGSF_EC'),
            variable = cms.string('pt')
        )
    ),
    mightGet = cms.optional.untracked.vstring,
    mvaOutput_normalization = cms.string(''),
    rawValues = cms.vstring(
        'discriminator',
        'category'
    ),
    toMultiplex = cms.InputTag("hpsPFTauDiscriminationByMVA6rawElectronRejection"),
    verbosity = cms.int32(0),
    workingPoints = cms.vstring(
        '_VLoose',
        '_Loose',
        '_Medium',
        '_Tight',
        '_VTight'
    )
)


process.hpsPFTauDiscriminationByMVA6rawElectronRejection = cms.EDProducer("PFRecoTauDiscriminationAgainstElectronMVA6",
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        decayMode = cms.PSet(
            Producer = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs"),
            cut = cms.double(0.5)
        )
    ),
    hgcalElectronIDs = cms.VInputTag(),
    isPhase2 = cms.bool(False),
    loadMVAfromDB = cms.bool(True),
    method = cms.string('BDTG'),
    mightGet = cms.optional.untracked.vstring,
    minMVANoEleMatchWOgWOgsfBL = cms.double(0),
    minMVANoEleMatchWOgWOgsfEC = cms.double(0),
    minMVANoEleMatchWgWOgsfBL = cms.double(0),
    minMVANoEleMatchWgWOgsfEC = cms.double(0),
    minMVAWOgWgsfBL = cms.double(0),
    minMVAWOgWgsfEC = cms.double(0),
    minMVAWgWgsfBL = cms.double(0),
    minMVAWgWgsfEC = cms.double(0),
    mvaName_NoEleMatch_wGwoGSF_BL = cms.string('RecoTauTag_antiElectronMVA_NoEleMatch_wGwoGSF_BL'),
    mvaName_NoEleMatch_wGwoGSF_EC = cms.string('RecoTauTag_antiElectronMVA_NoEleMatch_wGwoGSF_EC'),
    mvaName_NoEleMatch_woGwoGSF_BL = cms.string('RecoTauTag_antiElectronMVA_NoEleMatch_woGwoGSF_BL'),
    mvaName_NoEleMatch_woGwoGSF_EC = cms.string('RecoTauTag_antiElectronMVA_NoEleMatch_woGwoGSF_EC'),
    mvaName_wGwGSF_BL = cms.string('RecoTauTag_antiElectronMVA_wGwGSF_BL'),
    mvaName_wGwGSF_EC = cms.string('RecoTauTag_antiElectronMVA_wGwGSF_EC'),
    mvaName_woGwGSF_BL = cms.string('RecoTauTag_antiElectronMVA_woGwGSF_BL'),
    mvaName_woGwGSF_EC = cms.string('RecoTauTag_antiElectronMVA_woGwGSF_EC'),
    returnMVA = cms.bool(True),
    srcElectrons = cms.InputTag("gedGsfElectrons"),
    usePhiAtEcalEntranceExtrapolation = cms.bool(False),
    verbosity = cms.int32(0),
    vetoEcalCracks = cms.bool(False)
)


process.hpsPFTauDiscriminationByMediumElectronRejection = cms.EDProducer("PFRecoTauDiscriminationAgainstElectron",
    ApplyCut_BremCombined = cms.bool(False),
    ApplyCut_BremsRecoveryEOverPLead = cms.bool(False),
    ApplyCut_EOverPLead = cms.bool(False),
    ApplyCut_EcalCrackCut = cms.bool(True),
    ApplyCut_ElectronPreID = cms.bool(False),
    ApplyCut_ElectronPreID_2D = cms.bool(False),
    ApplyCut_EmFraction = cms.bool(False),
    ApplyCut_Hcal3x3OverPLead = cms.bool(False),
    ApplyCut_HcalMaxOverPLead = cms.bool(False),
    ApplyCut_HcalTotOverPLead = cms.bool(False),
    ApplyCut_PFElectronMVA = cms.bool(True),
    BremCombined_Fraction = cms.double(0.99),
    BremCombined_HOP = cms.double(0.1),
    BremCombined_Mass = cms.double(0.55),
    BremCombined_StripSize = cms.double(0.03),
    BremsRecoveryEOverPLead_maxValue = cms.double(1.8),
    BremsRecoveryEOverPLead_minValue = cms.double(0.8),
    EOverPLead_maxValue = cms.double(1.8),
    EOverPLead_minValue = cms.double(0.8),
    ElecPreID0_EOverPLead_maxValue = cms.double(0.95),
    ElecPreID0_HOverPLead_minValue = cms.double(0.05),
    ElecPreID1_EOverPLead_maxValue = cms.double(0.8),
    ElecPreID1_HOverPLead_minValue = cms.double(0.15),
    EmFraction_maxValue = cms.double(0.9),
    Hcal3x3OverPLead_minValue = cms.double(0.1),
    HcalMaxOverPLead_minValue = cms.double(0.1),
    HcalTotOverPLead_minValue = cms.double(0.1),
    PFElectronMVA_maxValue = cms.double(-0.1),
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and')
    )
)


process.hpsPFTauDiscriminationByMuonRejection3 = cms.EDProducer("PFRecoTauDiscriminationAgainstMuon2Container",
    IDWPdefinitions = cms.VPSet(
        cms.PSet(
            HoPMin = cms.double(0.2),
            IDname = cms.string('ByLooseMuonRejection3'),
            discriminatorOption = cms.string('custom'),
            doCaloMuonVeto = cms.bool(True),
            maxNumberOfHitsLast2Stations = cms.int32(-1),
            maxNumberOfMatches = cms.int32(1)
        ),
        cms.PSet(
            HoPMin = cms.double(0.2),
            IDname = cms.string('ByTightMuonRejection3'),
            discriminatorOption = cms.string('custom'),
            doCaloMuonVeto = cms.bool(True),
            maxNumberOfHitsLast2Stations = cms.int32(0),
            maxNumberOfMatches = cms.int32(1)
        )
    ),
    IDdefinitions = cms.VPSet(),
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and')
    ),
    dRmuonMatch = cms.double(0.3),
    dRmuonMatchLimitedToJetArea = cms.bool(False),
    maskHitsCSC = cms.vint32(0, 0, 0, 0),
    maskHitsDT = cms.vint32(0, 0, 0, 0),
    maskHitsRPC = cms.vint32(0, 0, 0, 0),
    maskMatchesCSC = cms.vint32(1, 0, 0, 0),
    maskMatchesDT = cms.vint32(0, 0, 0, 0),
    maskMatchesRPC = cms.vint32(0, 0, 0, 0),
    mightGet = cms.optional.untracked.vstring,
    minPtMatchedMuon = cms.double(5),
    srcMuons = cms.InputTag("muons"),
    verbosity = cms.int32(0)
)


process.hpsPFTauDiscriminationByTightElectronRejection = cms.EDProducer("PFRecoTauDiscriminationAgainstElectron",
    ApplyCut_BremCombined = cms.bool(True),
    ApplyCut_BremsRecoveryEOverPLead = cms.bool(False),
    ApplyCut_EOverPLead = cms.bool(False),
    ApplyCut_EcalCrackCut = cms.bool(True),
    ApplyCut_ElectronPreID = cms.bool(False),
    ApplyCut_ElectronPreID_2D = cms.bool(False),
    ApplyCut_EmFraction = cms.bool(False),
    ApplyCut_Hcal3x3OverPLead = cms.bool(False),
    ApplyCut_HcalMaxOverPLead = cms.bool(False),
    ApplyCut_HcalTotOverPLead = cms.bool(False),
    ApplyCut_PFElectronMVA = cms.bool(True),
    BremCombined_Fraction = cms.double(0.99),
    BremCombined_HOP = cms.double(0.1),
    BremCombined_Mass = cms.double(0.55),
    BremCombined_StripSize = cms.double(0.03),
    BremsRecoveryEOverPLead_maxValue = cms.double(1.8),
    BremsRecoveryEOverPLead_minValue = cms.double(0.8),
    EOverPLead_maxValue = cms.double(1.8),
    EOverPLead_minValue = cms.double(0.8),
    ElecPreID0_EOverPLead_maxValue = cms.double(0.95),
    ElecPreID0_HOverPLead_minValue = cms.double(0.05),
    ElecPreID1_EOverPLead_maxValue = cms.double(0.8),
    ElecPreID1_HOverPLead_minValue = cms.double(0.15),
    EmFraction_maxValue = cms.double(0.9),
    Hcal3x3OverPLead_minValue = cms.double(0.1),
    HcalMaxOverPLead_minValue = cms.double(0.1),
    HcalTotOverPLead_minValue = cms.double(0.1),
    PFElectronMVA_maxValue = cms.double(-0.1),
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and')
    )
)


process.hpsPFTauPrimaryVertexProducer = cms.EDProducer("PFTauPrimaryVertexProducer",
    Algorithm = cms.int32(0),
    ElectronTag = cms.InputTag(""),
    MuonTag = cms.InputTag(""),
    PFTauTag = cms.InputTag("hpsPFTauProducer"),
    PVTag = cms.InputTag("offlinePrimaryVertices"),
    RemoveElectronTracks = cms.bool(False),
    RemoveMuonTracks = cms.bool(False),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    cut = cms.string('pt > 18.0 & abs(eta) < 2.4'),
    discriminators = cms.VPSet(cms.PSet(
        discriminator = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs"),
        selectionCut = cms.double(0.5)
    )),
    mightGet = cms.optional.untracked.vstring,
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.2),
            maxDeltaZToLeadTrack = cms.double(-1.0),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.03),
            minGammaEt = cms.double(1.5),
            minTrackHits = cms.uint32(8),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(1.0),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        leadingTrkOrPFCandOption = cms.string('leadPFCand'),
        primaryVertexSrc = cms.InputTag("offlinePrimaryVertices"),
        pvFindingAlgo = cms.string('closestInDeltaZ'),
        recoverLeadingTrk = cms.bool(False),
        signalQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.4),
            maxDeltaZToLeadTrack = cms.double(-1.0),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(1.0),
            minNeutralHadronEt = cms.double(30.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        vertexTrackFiltering = cms.bool(False),
        vxAssocQualityCuts = cms.PSet(
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(1.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        )
    ),
    useBeamSpot = cms.bool(True),
    useSelectedTaus = cms.bool(False)
)


process.hpsPFTauProducer = cms.EDProducer("RecoTauPiZeroUnembedder",
    src = cms.InputTag("hpsPFTauProducerSansRefs")
)


process.hpsPFTauProducerSansRefs = cms.EDProducer("RecoTauCleaner",
    cleaners = cms.VPSet(
        cms.PSet(
            name = cms.string('Charge'),
            nprongs = cms.vuint32(1, 3),
            passForCharge = cms.int32(1),
            plugin = cms.string('RecoTauChargeCleanerPlugin'),
            selectionFailValue = cms.double(0),
            tolerance = cms.double(0)
        ),
        cms.PSet(
            name = cms.string('HPS_Select'),
            plugin = cms.string('RecoTauDiscriminantCleanerPlugin'),
            src = cms.InputTag("hpsSelectionDiscriminator"),
            tolerance = cms.double(0)
        ),
        cms.PSet(
            minTrackPt = cms.double(5.0),
            name = cms.string('killSoftTwoProngTaus'),
            plugin = cms.string('RecoTauSoftTwoProngTausCleanerPlugin'),
            tolerance = cms.double(0)
        ),
        cms.PSet(
            name = cms.string('ChargedHadronMultiplicity'),
            plugin = cms.string('RecoTauChargedHadronMultiplicityCleanerPlugin'),
            tolerance = cms.double(0)
        ),
        cms.PSet(
            name = cms.string('Pt'),
            plugin = cms.string('RecoTauStringCleanerPlugin'),
            selection = cms.string('leadCand().isNonnull()'),
            selectionFailValue = cms.double(1000.0),
            selectionPassFunction = cms.string('-pt()'),
            tolerance = cms.double(0.01)
        ),
        cms.PSet(
            name = cms.string('StripMultiplicity'),
            plugin = cms.string('RecoTauStringCleanerPlugin'),
            selection = cms.string('leadCand().isNonnull()'),
            selectionFailValue = cms.double(1000.0),
            selectionPassFunction = cms.string('-signalPiZeroCandidates().size()'),
            tolerance = cms.double(0)
        ),
        cms.PSet(
            name = cms.string('CombinedIsolation'),
            plugin = cms.string('RecoTauStringCleanerPlugin'),
            selection = cms.string('leadCand().isNonnull()'),
            selectionFailValue = cms.double(1000.0),
            selectionPassFunction = cms.string('isolationPFChargedHadrCandsPtSum() + isolationPFGammaCandsEtSum()'),
            tolerance = cms.double(0)
        )
    ),
    outputSelection = cms.string(''),
    src = cms.InputTag("combinatoricRecoTaus"),
    verbosity = cms.int32(0)
)


process.hpsPFTauSecondaryVertexProducer = cms.EDProducer("PFTauSecondaryVertexProducer",
    PFTauTag = cms.InputTag("hpsPFTauProducer"),
    mightGet = cms.optional.untracked.vstring
)


process.hpsPFTauTransverseImpactParameters = cms.EDProducer("PFTauTransverseImpactParameters",
    PFTauPVATag = cms.InputTag("hpsPFTauPrimaryVertexProducer"),
    PFTauSVATag = cms.InputTag("hpsPFTauSecondaryVertexProducer"),
    PFTauTag = cms.InputTag("hpsPFTauProducer"),
    mightGet = cms.optional.untracked.vstring,
    useFullCalculation = cms.bool(True)
)


process.hpsSelectionDiscriminator = cms.EDProducer("PFRecoTauDiscriminationByHPSSelection",
    PFTauProducer = cms.InputTag("combinatoricRecoTaus"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and')
    ),
    decayModes = cms.VPSet(
        cms.PSet(
            applyBendCorrection = cms.PSet(
                eta = cms.bool(True),
                mass = cms.bool(True),
                phi = cms.bool(True)
            ),
            assumeStripMass = cms.double(-1.0),
            maxMass = cms.string('1.'),
            maxPi0Mass = cms.double(1000000000.0),
            minMass = cms.double(-1000.0),
            minPi0Mass = cms.double(-1000.0),
            nCharged = cms.uint32(1),
            nChargedPFCandsMin = cms.uint32(1),
            nPiZeros = cms.uint32(0),
            nTracksMin = cms.uint32(1)
        ),
        cms.PSet(
            applyBendCorrection = cms.PSet(
                eta = cms.bool(True),
                mass = cms.bool(True),
                phi = cms.bool(True)
            ),
            assumeStripMass = cms.double(0.1349),
            maxMass = cms.string('max(1.3, min(1.3*sqrt(pt/100.), 4.2))'),
            maxPi0Mass = cms.double(1000000000.0),
            minMass = cms.double(0.3),
            minPi0Mass = cms.double(-1000.0),
            nCharged = cms.uint32(1),
            nChargedPFCandsMin = cms.uint32(1),
            nPiZeros = cms.uint32(1),
            nTracksMin = cms.uint32(1)
        ),
        cms.PSet(
            applyBendCorrection = cms.PSet(
                eta = cms.bool(True),
                mass = cms.bool(True),
                phi = cms.bool(True)
            ),
            assumeStripMass = cms.double(0.0),
            maxMass = cms.string('max(1.2, min(1.2*sqrt(pt/100.), 4.0))'),
            maxPi0Mass = cms.double(0.2),
            minMass = cms.double(0.4),
            minPi0Mass = cms.double(0.05),
            nCharged = cms.uint32(1),
            nChargedPFCandsMin = cms.uint32(1),
            nPiZeros = cms.uint32(2),
            nTracksMin = cms.uint32(1)
        ),
        cms.PSet(
            applyBendCorrection = cms.PSet(
                eta = cms.bool(False),
                mass = cms.bool(False),
                phi = cms.bool(False)
            ),
            assumeStripMass = cms.double(-1.0),
            maxMass = cms.string('1.2'),
            maxPi0Mass = cms.double(1000000000.0),
            minMass = cms.double(0.0),
            minPi0Mass = cms.double(-1000.0),
            nCharged = cms.uint32(2),
            nChargedPFCandsMin = cms.uint32(1),
            nPiZeros = cms.uint32(0),
            nTracksMin = cms.uint32(2)
        ),
        cms.PSet(
            applyBendCorrection = cms.PSet(
                eta = cms.bool(False),
                mass = cms.bool(False),
                phi = cms.bool(False)
            ),
            assumeStripMass = cms.double(-1.0),
            maxMass = cms.string('max(1.2, min(1.2*sqrt(pt/100.), 4.0))'),
            maxPi0Mass = cms.double(1000000000.0),
            minMass = cms.double(0.0),
            minPi0Mass = cms.double(-1000.0),
            nCharged = cms.uint32(2),
            nChargedPFCandsMin = cms.uint32(1),
            nPiZeros = cms.uint32(1),
            nTracksMin = cms.uint32(2)
        ),
        cms.PSet(
            applyBendCorrection = cms.PSet(
                eta = cms.bool(False),
                mass = cms.bool(False),
                phi = cms.bool(False)
            ),
            assumeStripMass = cms.double(-1.0),
            maxMass = cms.string('1.5'),
            maxPi0Mass = cms.double(1000000000.0),
            minMass = cms.double(0.8),
            minPi0Mass = cms.double(-1000.0),
            nCharged = cms.uint32(3),
            nChargedPFCandsMin = cms.uint32(1),
            nPiZeros = cms.uint32(0),
            nTracksMin = cms.uint32(2)
        ),
        cms.PSet(
            applyBendCorrection = cms.PSet(
                eta = cms.bool(False),
                mass = cms.bool(False),
                phi = cms.bool(False)
            ),
            assumeStripMass = cms.double(-1.0),
            maxMass = cms.string('1.6'),
            maxPi0Mass = cms.double(1000000000.0),
            minMass = cms.double(0.9),
            minPi0Mass = cms.double(-1000.0),
            nCharged = cms.uint32(3),
            nChargedPFCandsMin = cms.uint32(1),
            nPiZeros = cms.uint32(1),
            nTracksMin = cms.uint32(2)
        )
    ),
    matchingCone = cms.double(0.5),
    minPixelHits = cms.int32(1),
    minTauPt = cms.double(0.0),
    requireTauChargedHadronsToBeChargedPFCands = cms.bool(False),
    verbosity = cms.int32(0)
)


process.initialStep = cms.EDProducer("TrackTfClassifier",
    beamspot = cms.InputTag("offlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mightGet = cms.optional.untracked.vstring,
    mva = cms.PSet(
        batchSize = cms.int32(16),
        tfDnnLabel = cms.string('trackSelectionTf')
    ),
    qualityCuts = cms.vdouble(-0.35, 0.1, 0.28),
    src = cms.InputTag("initialStepTracks"),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.initialStepClassifier1 = cms.EDProducer("TrackMVAClassifierPrompt",
    beamspot = cms.InputTag("offlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mightGet = cms.optional.untracked.vstring,
    mva = cms.PSet(
        GBRForestFileName = cms.string(''),
        GBRForestLabel = cms.string('MVASelectorIter0_13TeV')
    ),
    qualityCuts = cms.vdouble(-0.9, -0.8, -0.7),
    src = cms.InputTag("initialStepTracks"),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.initialStepClassifier2 = cms.EDProducer("TrackMVAClassifierDetached",
    beamspot = cms.InputTag("offlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mightGet = cms.optional.untracked.vstring,
    mva = cms.PSet(
        GBRForestFileName = cms.string(''),
        GBRForestLabel = cms.string('MVASelectorIter3_13TeV')
    ),
    qualityCuts = cms.vdouble(-0.5, 0.0, 0.5),
    src = cms.InputTag("initialStepTracks"),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.initialStepClassifier3 = cms.EDProducer("TrackTfClassifier",
    beamspot = cms.InputTag("offlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mightGet = cms.optional.untracked.vstring,
    mva = cms.PSet(
        batchSize = cms.int32(16),
        tfDnnLabel = cms.string('trackSelectionTf_CKF')
    ),
    qualityCuts = cms.vdouble(-0.23, 0.15, 0.41),
    src = cms.InputTag("initialStepTracks"),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.initialStepHitDoublets = cms.EDProducer("HitPairEDProducer",
    clusterCheck = cms.InputTag("trackerClusterCheck"),
    layerPairs = cms.vuint32(0, 1, 2),
    maxElement = cms.uint32(50000000),
    maxElementTotal = cms.uint32(50000000),
    mightGet = cms.optional.untracked.vstring,
    produceIntermediateHitDoublets = cms.bool(True),
    produceSeedingHitSets = cms.bool(False),
    seedingLayers = cms.InputTag("initialStepSeedLayers"),
    trackingRegions = cms.InputTag("initialStepTrackingRegions"),
    trackingRegionsSeedingLayers = cms.InputTag("")
)


process.initialStepHitDoubletsPreSplitting = cms.EDProducer("HitPairEDProducer",
    clusterCheck = cms.InputTag("trackerClusterCheckPreSplitting"),
    layerPairs = cms.vuint32(0, 1, 2),
    maxElement = cms.uint32(50000000),
    maxElementTotal = cms.uint32(50000000),
    mightGet = cms.optional.untracked.vstring,
    produceIntermediateHitDoublets = cms.bool(True),
    produceSeedingHitSets = cms.bool(False),
    seedingLayers = cms.InputTag("initialStepSeedLayersPreSplitting"),
    trackingRegions = cms.InputTag("initialStepTrackingRegionsPreSplitting"),
    trackingRegionsSeedingLayers = cms.InputTag("")
)


process.initialStepHitQuadruplets = cms.EDProducer("CAHitQuadrupletEDProducer",
    CAHardPtCut = cms.double(0),
    CAOnlyOneLastHitPerLayerFilter = cms.optional.bool,
    CAPhiCut = cms.double(0.2),
    CAPhiCut_byTriplets = cms.VPSet(cms.PSet(
        cut = cms.double(-1),
        seedingLayers = cms.string('')
    )),
    CAThetaCut = cms.double(0.0012),
    CAThetaCut_byTriplets = cms.VPSet(cms.PSet(
        cut = cms.double(-1),
        seedingLayers = cms.string('')
    )),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('LowPtClusterShapeSeedComparitor'),
        clusterShapeCacheSrc = cms.InputTag("siPixelClusterShapeCache"),
        clusterShapeHitFilter = cms.string('ClusterShapeHitFilter')
    ),
    doublets = cms.InputTag("initialStepHitDoublets"),
    extraHitRPhitolerance = cms.double(0.032),
    fitFastCircle = cms.bool(True),
    fitFastCircleChi2Cut = cms.bool(True),
    maxChi2 = cms.PSet(
        enabled = cms.bool(True),
        pt1 = cms.double(0.7),
        pt2 = cms.double(2),
        value1 = cms.double(200),
        value2 = cms.double(50)
    ),
    mightGet = cms.optional.untracked.vstring,
    useBendingCorrection = cms.bool(True)
)


process.initialStepHitQuadrupletsPreSplitting = cms.EDProducer("CAHitQuadrupletEDProducer",
    CAHardPtCut = cms.double(0),
    CAOnlyOneLastHitPerLayerFilter = cms.optional.bool,
    CAPhiCut = cms.double(0.2),
    CAPhiCut_byTriplets = cms.VPSet(cms.PSet(
        cut = cms.double(-1),
        seedingLayers = cms.string('')
    )),
    CAThetaCut = cms.double(0.0012),
    CAThetaCut_byTriplets = cms.VPSet(cms.PSet(
        cut = cms.double(-1),
        seedingLayers = cms.string('')
    )),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('LowPtClusterShapeSeedComparitor'),
        clusterShapeCacheSrc = cms.InputTag("siPixelClusterShapeCachePreSplitting"),
        clusterShapeHitFilter = cms.string('ClusterShapeHitFilter')
    ),
    doublets = cms.InputTag("initialStepHitDoubletsPreSplitting"),
    extraHitRPhitolerance = cms.double(0.032),
    fitFastCircle = cms.bool(True),
    fitFastCircleChi2Cut = cms.bool(True),
    maxChi2 = cms.PSet(
        enabled = cms.bool(True),
        pt1 = cms.double(0.7),
        pt2 = cms.double(2),
        value1 = cms.double(200),
        value2 = cms.double(50)
    ),
    mightGet = cms.optional.untracked.vstring,
    useBendingCorrection = cms.bool(True)
)


process.initialStepHitTriplets = cms.EDProducer("PixelTripletHLTEDProducer",
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('LowPtClusterShapeSeedComparitor'),
        clusterShapeCacheSrc = cms.InputTag("siPixelClusterShapeCache"),
        clusterShapeHitFilter = cms.string('ClusterShapeHitFilter')
    ),
    doublets = cms.InputTag("initialStepHitDoublets"),
    extraHitRPhitolerance = cms.double(0.032),
    extraHitRZtolerance = cms.double(0.037),
    maxElement = cms.uint32(1000000),
    mightGet = cms.optional.untracked.vstring,
    phiPreFiltering = cms.double(0.3),
    produceIntermediateHitTriplets = cms.bool(False),
    produceSeedingHitSets = cms.bool(True),
    useBending = cms.bool(True),
    useFixedPreFiltering = cms.bool(False),
    useMultScattering = cms.bool(True)
)


process.initialStepHitTripletsPreSplitting = cms.EDProducer("PixelTripletHLTEDProducer",
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('LowPtClusterShapeSeedComparitor'),
        clusterShapeCacheSrc = cms.InputTag("siPixelClusterShapeCachePreSplitting"),
        clusterShapeHitFilter = cms.string('ClusterShapeHitFilter')
    ),
    doublets = cms.InputTag("initialStepHitDoubletsPreSplitting"),
    extraHitRPhitolerance = cms.double(0.032),
    extraHitRZtolerance = cms.double(0.037),
    maxElement = cms.uint32(1000000),
    mightGet = cms.optional.untracked.vstring,
    phiPreFiltering = cms.double(0.3),
    produceIntermediateHitTriplets = cms.bool(False),
    produceSeedingHitSets = cms.bool(True),
    useBending = cms.bool(True),
    useFixedPreFiltering = cms.bool(False),
    useMultScattering = cms.bool(True)
)


process.initialStepSeedClusterMask = cms.EDProducer("SeedClusterRemover",
    Common = cms.PSet(
        maxChi2 = cms.double(9.0)
    ),
    oldClusterRemovalInfo = cms.InputTag("pixelLessStepClusters"),
    pixelClusters = cms.InputTag("siPixelClusters"),
    stripClusters = cms.InputTag("siStripClusters"),
    trajectories = cms.InputTag("initialStepSeeds")
)


process.initialStepSeedLayers = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(
        HitProducer = cms.string('siPixelRecHits'),
        TTRHBuilder = cms.string('WithTrackAngle')
    ),
    FPix = cms.PSet(
        HitProducer = cms.string('siPixelRecHits'),
        TTRHBuilder = cms.string('WithTrackAngle')
    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(

    ),
    TIB = cms.PSet(

    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(

    ),
    layerList = cms.vstring(
        'BPix1+BPix2+BPix3+BPix4',
        'BPix1+BPix2+BPix3+FPix1_pos',
        'BPix1+BPix2+BPix3+FPix1_neg',
        'BPix1+BPix2+FPix1_pos+FPix2_pos',
        'BPix1+BPix2+FPix1_neg+FPix2_neg',
        'BPix1+FPix1_pos+FPix2_pos+FPix3_pos',
        'BPix1+FPix1_neg+FPix2_neg+FPix3_neg'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.initialStepSeedLayersPreSplitting = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(
        HitProducer = cms.string('siPixelRecHitsPreSplitting'),
        TTRHBuilder = cms.string('WithTrackAngle')
    ),
    FPix = cms.PSet(
        HitProducer = cms.string('siPixelRecHitsPreSplitting'),
        TTRHBuilder = cms.string('WithTrackAngle')
    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(

    ),
    TIB = cms.PSet(

    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(

    ),
    layerList = cms.vstring(
        'BPix1+BPix2+BPix3+BPix4',
        'BPix1+BPix2+BPix3+FPix1_pos',
        'BPix1+BPix2+BPix3+FPix1_neg',
        'BPix1+BPix2+FPix1_pos+FPix2_pos',
        'BPix1+BPix2+FPix1_neg+FPix2_neg',
        'BPix1+FPix1_pos+FPix2_pos+FPix3_pos',
        'BPix1+FPix1_neg+FPix2_neg+FPix3_neg'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.initialStepSeeds = cms.EDProducer("SeedCreatorFromRegionConsecutiveHitsTripletOnlyEDProducer",
    MinOneOverPtError = cms.double(1),
    OriginTransverseErrorMultiplier = cms.double(1),
    SeedComparitorPSet = cms.PSet(
        ClusterShapeCacheSrc = cms.InputTag("siPixelClusterShapeCache"),
        ClusterShapeHitFilterName = cms.string('ClusterShapeHitFilter'),
        ComponentName = cms.string('PixelClusterShapeSeedComparitor'),
        FilterAtHelixStage = cms.bool(False),
        FilterPixelHits = cms.bool(True),
        FilterStripHits = cms.bool(False)
    ),
    SeedMomentumForBOFF = cms.double(5),
    TTRHBuilder = cms.string('WithTrackAngle'),
    forceKinematicWithRegionDirection = cms.bool(False),
    magneticField = cms.string('ParabolicMf'),
    mightGet = cms.optional.untracked.vstring,
    propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    seedingHitSets = cms.InputTag("initialStepHitQuadruplets")
)


process.initialStepSeedsPreSplitting = cms.EDProducer("SeedCreatorFromRegionConsecutiveHitsEDProducer",
    MinOneOverPtError = cms.double(1),
    OriginTransverseErrorMultiplier = cms.double(1),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('none')
    ),
    SeedMomentumForBOFF = cms.double(5),
    TTRHBuilder = cms.string('WithTrackAngle'),
    forceKinematicWithRegionDirection = cms.bool(False),
    magneticField = cms.string('ParabolicMf'),
    mightGet = cms.optional.untracked.vstring,
    propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    seedingHitSets = cms.InputTag("initialStepHitQuadrupletsPreSplitting")
)


process.initialStepSelector = cms.EDProducer("MultiTrackSelector",
    GBRForestLabel = cms.string('MVASelectorIter0'),
    beamspot = cms.InputTag("offlineBeamSpot"),
    src = cms.InputTag("initialStepTracks"),
    trackSelectors = cms.VPSet(
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(True),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(1.6),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(0.55, 4.0),
            d0_par2 = cms.vdouble(0.55, 4.0),
            dz_par1 = cms.vdouble(0.65, 4.0),
            dz_par2 = cms.vdouble(0.45, 4.0),
            keepAllTracks = cms.bool(False),
            maxNumberLostLayers = cms.uint32(999),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(0),
            minNumberLayers = cms.uint32(0),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(0),
            nSigmaZ = cms.double(4.0),
            name = cms.string('initialStepLoose'),
            preFilterName = cms.string(''),
            qualityBit = cms.string('loose'),
            res_par = cms.vdouble(0.003, 0.01),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        ),
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(True),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(0.7),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(0.3, 4.0),
            d0_par2 = cms.vdouble(0.4, 4.0),
            dz_par1 = cms.vdouble(0.35, 4.0),
            dz_par2 = cms.vdouble(0.4, 4.0),
            keepAllTracks = cms.bool(True),
            maxNumberLostLayers = cms.uint32(2),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(3),
            minNumberLayers = cms.uint32(3),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(0),
            nSigmaZ = cms.double(4.0),
            name = cms.string('initialStepTight'),
            preFilterName = cms.string('initialStepLoose'),
            qualityBit = cms.string('tight'),
            res_par = cms.vdouble(0.003, 0.01),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        ),
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(True),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(0.7),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(0.3, 4.0),
            d0_par2 = cms.vdouble(0.4, 4.0),
            dz_par1 = cms.vdouble(0.35, 4.0),
            dz_par2 = cms.vdouble(0.4, 4.0),
            keepAllTracks = cms.bool(True),
            maxNumberLostLayers = cms.uint32(2),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(3),
            minNumberLayers = cms.uint32(3),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(0),
            nSigmaZ = cms.double(4.0),
            name = cms.string('QualityMasks'),
            preFilterName = cms.string('initialStepTight'),
            qualityBit = cms.string('highPurity'),
            res_par = cms.vdouble(0.003, 0.001),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        )
    ),
    useAnyMVA = cms.bool(False),
    useVertices = cms.bool(True),
    useVtxError = cms.bool(False),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.initialStepTrackCandidates = cms.EDProducer("MkFitOutputConverter",
    batchSize = cms.int32(16),
    candMVASel = cms.bool(False),
    candWP = cms.double(0),
    doErrorRescale = cms.bool(True),
    mightGet = cms.optional.untracked.vstring,
    mkFitEventOfHits = cms.InputTag("mkFitEventOfHits"),
    mkFitPixelHits = cms.InputTag("mkFitSiPixelHits"),
    mkFitSeeds = cms.InputTag("initialStepTrackCandidatesMkFitSeeds"),
    mkFitStripHits = cms.InputTag("mkFitSiStripHits"),
    propagatorAlong = cms.ESInputTag("","PropagatorWithMaterial"),
    propagatorOpposite = cms.ESInputTag("","PropagatorWithMaterialOpposite"),
    qualityMaxInvPt = cms.double(100),
    qualityMaxPosErr = cms.double(100),
    qualityMaxR = cms.double(120),
    qualityMaxZ = cms.double(280),
    qualityMinTheta = cms.double(0.01),
    qualitySignPt = cms.bool(True),
    seeds = cms.InputTag("initialStepSeeds"),
    tfDnnLabel = cms.string('trackSelectionTf'),
    tracks = cms.InputTag("initialStepTrackCandidatesMkFit"),
    ttrhBuilder = cms.ESInputTag("","WithTrackAngle")
)


process.initialStepTrackCandidatesMkFit = cms.EDProducer("MkFitProducer",
    backwardFitInCMSSW = cms.bool(False),
    buildingRoutine = cms.string('cloneEngine'),
    clustersToSkip = cms.InputTag(""),
    config = cms.ESInputTag("","initialStepTrackCandidatesMkFitConfig"),
    eventOfHits = cms.InputTag("mkFitEventOfHits"),
    limitConcurrency = cms.untracked.bool(False),
    mightGet = cms.optional.untracked.vstring,
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutLoose')
    ),
    mkFitSilent = cms.untracked.bool(True),
    pixelHits = cms.InputTag("mkFitSiPixelHits"),
    removeDuplicates = cms.bool(True),
    seedCleaning = cms.bool(True),
    seeds = cms.InputTag("initialStepTrackCandidatesMkFitSeeds"),
    stripHits = cms.InputTag("mkFitSiStripHits")
)


process.initialStepTrackCandidatesMkFitPreSplitting = cms.EDProducer("MkFitProducer",
    backwardFitInCMSSW = cms.bool(False),
    buildingRoutine = cms.string('cloneEngine'),
    clustersToSkip = cms.InputTag(""),
    config = cms.ESInputTag("","initialStepTrackCandidatesMkFitConfigPreSplitting"),
    eventOfHits = cms.InputTag("mkFitEventOfHitsPreSplitting"),
    limitConcurrency = cms.untracked.bool(False),
    mightGet = cms.optional.untracked.vstring,
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutLoose')
    ),
    mkFitSilent = cms.untracked.bool(True),
    pixelHits = cms.InputTag("mkFitSiPixelHitsPreSplitting"),
    removeDuplicates = cms.bool(True),
    seedCleaning = cms.bool(True),
    seeds = cms.InputTag("initialStepTrackCandidatesMkFitSeedsPreSplitting"),
    stripHits = cms.InputTag("mkFitSiStripHits")
)


process.initialStepTrackCandidatesMkFitSeeds = cms.EDProducer("MkFitSeedConverter",
    maxNSeeds = cms.uint32(500000),
    mightGet = cms.optional.untracked.vstring,
    seeds = cms.InputTag("initialStepSeeds"),
    ttrhBuilder = cms.ESInputTag("","WithTrackAngle")
)


process.initialStepTrackCandidatesMkFitSeedsPreSplitting = cms.EDProducer("MkFitSeedConverter",
    maxNSeeds = cms.uint32(500000),
    mightGet = cms.optional.untracked.vstring,
    seeds = cms.InputTag("initialStepSeedsPreSplitting"),
    ttrhBuilder = cms.ESInputTag("","WithTrackAngle")
)


process.initialStepTrackCandidatesPreSplitting = cms.EDProducer("MkFitOutputConverter",
    batchSize = cms.int32(16),
    candMVASel = cms.bool(False),
    candWP = cms.double(0),
    doErrorRescale = cms.bool(True),
    mightGet = cms.optional.untracked.vstring,
    mkFitEventOfHits = cms.InputTag("mkFitEventOfHitsPreSplitting"),
    mkFitPixelHits = cms.InputTag("mkFitSiPixelHitsPreSplitting"),
    mkFitSeeds = cms.InputTag("initialStepTrackCandidatesMkFitSeedsPreSplitting"),
    mkFitStripHits = cms.InputTag("mkFitSiStripHits"),
    propagatorAlong = cms.ESInputTag("","PropagatorWithMaterial"),
    propagatorOpposite = cms.ESInputTag("","PropagatorWithMaterialOpposite"),
    qualityMaxInvPt = cms.double(100),
    qualityMaxPosErr = cms.double(100),
    qualityMaxR = cms.double(120),
    qualityMaxZ = cms.double(280),
    qualityMinTheta = cms.double(0.01),
    qualitySignPt = cms.bool(True),
    seeds = cms.InputTag("initialStepSeedsPreSplitting"),
    tfDnnLabel = cms.string('trackSelectionTf'),
    tracks = cms.InputTag("initialStepTrackCandidatesMkFitPreSplitting"),
    ttrhBuilder = cms.ESInputTag("","WithTrackAngle")
)


process.initialStepTrackRefsForJets = cms.EDProducer("ChargedRefCandidateProducer",
    particleType = cms.string('pi+'),
    src = cms.InputTag("initialStepTracks")
)


process.initialStepTrackRefsForJetsPreSplitting = cms.EDProducer("ChargedRefCandidateProducer",
    particleType = cms.string('pi+'),
    src = cms.InputTag("initialStepTracksPreSplitting")
)


process.initialStepTrackingRegions = cms.EDProducer("GlobalTrackingRegionFromBeamSpotEDProducer",
    RegionPSet = cms.PSet(
        beamSpot = cms.InputTag("offlineBeamSpot"),
        nSigmaZ = cms.double(4),
        originHalfLength = cms.double(0),
        originRadius = cms.double(0.02),
        precise = cms.bool(True),
        ptMin = cms.double(0.5),
        useMultipleScattering = cms.bool(False)
    ),
    mightGet = cms.optional.untracked.vstring
)


process.initialStepTrackingRegionsPreSplitting = cms.EDProducer("GlobalTrackingRegionFromBeamSpotEDProducer",
    RegionPSet = cms.PSet(
        beamSpot = cms.InputTag("offlineBeamSpot"),
        nSigmaZ = cms.double(4),
        originHalfLength = cms.double(0),
        originRadius = cms.double(0.02),
        precise = cms.bool(True),
        ptMin = cms.double(0.5),
        useMultipleScattering = cms.bool(False)
    ),
    mightGet = cms.optional.untracked.vstring
)


process.initialStepTracks = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('initialStep'),
    Fitter = cms.string('FlexibleKFFittingSmoother'),
    GeometricInnerState = cms.bool(False),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    Propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    TTRHBuilder = cms.string('WithAngleAndTemplateWithoutProbQ'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("initialStepTrackCandidates"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(True)
)


process.initialStepTracksPreSplitting = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('initialStep'),
    Fitter = cms.string('FlexibleKFFittingSmoother'),
    GeometricInnerState = cms.bool(False),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEventPreSplitting"),
    NavigationSchool = cms.string(''),
    Propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    TTRHBuilder = cms.string('WithAngleAndTemplateWithoutProbQ'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("initialStepTrackCandidatesPreSplitting"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(True)
)


process.jetCoreClusterSplitter = cms.EDProducer("JetCoreClusterSplitter",
    centralMIPCharge = cms.double(26000),
    chargeFractionMin = cms.double(2.0),
    chargePerUnit = cms.double(2000),
    cores = cms.InputTag("ak5CaloJets"),
    deltaRmax = cms.double(0.05),
    forceXError = cms.double(100),
    forceYError = cms.double(150),
    fractionalWidth = cms.double(0.4),
    pixelCPE = cms.string('PixelCPEGeneric'),
    pixelClusters = cms.InputTag("siPixelCluster"),
    ptMin = cms.double(200),
    verbose = cms.bool(False),
    vertices = cms.InputTag("offlinePrimaryVertices")
)


process.jetCoreRegionalStep = cms.EDProducer("TrackTfClassifier",
    beamspot = cms.InputTag("offlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mightGet = cms.optional.untracked.vstring,
    mva = cms.PSet(
        batchSize = cms.int32(16),
        tfDnnLabel = cms.string('trackSelectionTf')
    ),
    qualityCuts = cms.vdouble(-0.62, -0.49, 0.02),
    src = cms.InputTag("jetCoreRegionalStepTracks"),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.jetCoreRegionalStepBarrel = cms.EDProducer("TrackTfClassifier",
    beamspot = cms.InputTag("offlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mightGet = cms.optional.untracked.vstring,
    mva = cms.PSet(
        batchSize = cms.int32(16),
        tfDnnLabel = cms.string('trackSelectionTf')
    ),
    qualityCuts = cms.vdouble(-0.62, -0.49, 0.02),
    src = cms.InputTag("jetCoreRegionalStepBarrelTracks"),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.jetCoreRegionalStepBarrelTrackCandidates = cms.EDProducer("CkfTrackCandidateMaker",
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    RedundantSeedCleaner = cms.string('CachingSeedCleanerBySharedInput'),
    TrajectoryBuilderPSet = cms.PSet(
        refToPSet_ = cms.string('jetCoreRegionalStepBarrelTrajectoryBuilder')
    ),
    TrajectoryCleaner = cms.string('jetCoreRegionalStepDeepCoreTrajectoryCleaner'),
    TransientInitialStateEstimatorParameters = cms.PSet(
        numberMeasurementsForFit = cms.int32(4),
        propagatorAlongTISE = cms.string('PropagatorWithMaterial'),
        propagatorOppositeTISE = cms.string('PropagatorWithMaterialOpposite')
    ),
    cleanTrajectoryAfterInOut = cms.bool(True),
    clustersToSkip = cms.InputTag(""),
    doSeedingRegionRebuilding = cms.bool(True),
    maxNSeeds = cms.uint32(500000),
    maxSeedsBeforeCleaning = cms.uint32(10000),
    mightGet = cms.optional.untracked.vstring,
    numHitsForSeedCleaner = cms.int32(4),
    onlyPixelHitsForSeedCleaner = cms.bool(False),
    phase2clustersToSkip = cms.InputTag(""),
    reverseTrajectories = cms.bool(False),
    src = cms.InputTag("jetCoreRegionalStepSeedsBarrel"),
    useHitsSplitting = cms.bool(True)
)


process.jetCoreRegionalStepBarrelTracks = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('jetCoreRegionalStep'),
    Fitter = cms.string('FlexibleKFFittingSmoother'),
    GeometricInnerState = cms.bool(False),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    Propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    TTRHBuilder = cms.string('WithAngleAndTemplateWithoutProbQ'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("jetCoreRegionalStepBarrelTrackCandidates"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(True)
)


process.jetCoreRegionalStepEndcap = cms.EDProducer("TrackTfClassifier",
    beamspot = cms.InputTag("offlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mightGet = cms.optional.untracked.vstring,
    mva = cms.PSet(
        batchSize = cms.int32(16),
        tfDnnLabel = cms.string('trackSelectionTf')
    ),
    qualityCuts = cms.vdouble(-0.62, -0.49, 0.02),
    src = cms.InputTag("jetCoreRegionalStepEndcapTracks"),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.jetCoreRegionalStepEndcapHitDoublets = cms.EDProducer("HitPairEDProducer",
    clusterCheck = cms.InputTag("trackerClusterCheck"),
    layerPairs = cms.vuint32(0),
    maxElement = cms.uint32(1000000),
    maxElementTotal = cms.uint32(12000000),
    mightGet = cms.optional.untracked.vstring,
    produceIntermediateHitDoublets = cms.bool(False),
    produceSeedingHitSets = cms.bool(True),
    seedingLayers = cms.InputTag("jetCoreRegionalStepSeedLayers"),
    trackingRegions = cms.InputTag("jetCoreRegionalStepEndcapTrackingRegions"),
    trackingRegionsSeedingLayers = cms.InputTag("")
)


process.jetCoreRegionalStepEndcapTrackCandidates = cms.EDProducer("CkfTrackCandidateMaker",
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    RedundantSeedCleaner = cms.string('CachingSeedCleanerBySharedInput'),
    TrajectoryBuilderPSet = cms.PSet(
        refToPSet_ = cms.string('jetCoreRegionalStepEndcapTrajectoryBuilder')
    ),
    TrajectoryCleaner = cms.string('TrajectoryCleanerBySharedHits'),
    TransientInitialStateEstimatorParameters = cms.PSet(
        numberMeasurementsForFit = cms.int32(4),
        propagatorAlongTISE = cms.string('PropagatorWithMaterial'),
        propagatorOppositeTISE = cms.string('PropagatorWithMaterialOpposite')
    ),
    cleanTrajectoryAfterInOut = cms.bool(True),
    clustersToSkip = cms.InputTag(""),
    doSeedingRegionRebuilding = cms.bool(True),
    maxNSeeds = cms.uint32(500000),
    maxSeedsBeforeCleaning = cms.uint32(10000),
    mightGet = cms.optional.untracked.vstring,
    numHitsForSeedCleaner = cms.int32(4),
    onlyPixelHitsForSeedCleaner = cms.bool(False),
    phase2clustersToSkip = cms.InputTag(""),
    reverseTrajectories = cms.bool(False),
    src = cms.InputTag("jetCoreRegionalStepSeedsEndcap"),
    useHitsSplitting = cms.bool(True)
)


process.jetCoreRegionalStepEndcapTrackingRegions = cms.EDProducer("TauRegionalPixelSeedTrackingRegionEDProducer",
    RegionPSet = cms.PSet(
        JetSrc = cms.InputTag("jetsForCoreTrackingEndcap"),
        deltaEtaRegion = cms.double(0.2),
        deltaPhiRegion = cms.double(0.2),
        howToUseMeasurementTracker = cms.string('Never'),
        measurementTrackerName = cms.InputTag("MeasurementTrackerEvent"),
        originHalfLength = cms.double(0.2),
        originRadius = cms.double(0.2),
        ptMin = cms.double(10),
        searchOpt = cms.bool(False),
        vertexSrc = cms.InputTag("firstStepGoodPrimaryVertices")
    ),
    mightGet = cms.optional.untracked.vstring
)


process.jetCoreRegionalStepEndcapTracks = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('jetCoreRegionalStep'),
    Fitter = cms.string('FlexibleKFFittingSmoother'),
    GeometricInnerState = cms.bool(False),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    Propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    TTRHBuilder = cms.string('WithAngleAndTemplateWithoutProbQ'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("jetCoreRegionalStepEndcapTrackCandidates"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(True)
)


process.jetCoreRegionalStepHitDoublets = cms.EDProducer("HitPairEDProducer",
    clusterCheck = cms.InputTag("trackerClusterCheck"),
    layerPairs = cms.vuint32(0),
    maxElement = cms.uint32(1000000),
    maxElementTotal = cms.uint32(12000000),
    mightGet = cms.optional.untracked.vstring,
    produceIntermediateHitDoublets = cms.bool(False),
    produceSeedingHitSets = cms.bool(True),
    seedingLayers = cms.InputTag("jetCoreRegionalStepSeedLayers"),
    trackingRegions = cms.InputTag("jetCoreRegionalStepTrackingRegions"),
    trackingRegionsSeedingLayers = cms.InputTag("")
)


process.jetCoreRegionalStepSeedLayers = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(
        HitProducer = cms.string('siPixelRecHits'),
        TTRHBuilder = cms.string('WithTrackAngle'),
        hitErrorRPhi = cms.double(0.0027),
        hitErrorRZ = cms.double(0.006),
        useErrorsFromParam = cms.bool(True)
    ),
    FPix = cms.PSet(
        HitProducer = cms.string('siPixelRecHits'),
        TTRHBuilder = cms.string('WithTrackAngle'),
        hitErrorRPhi = cms.double(0.0051),
        hitErrorRZ = cms.double(0.0036),
        useErrorsFromParam = cms.bool(True)
    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(

    ),
    TIB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit")
    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(

    ),
    layerList = cms.vstring(
        'BPix1+BPix2',
        'BPix1+BPix3',
        'BPix1+BPix4',
        'BPix2+BPix3',
        'BPix2+BPix4',
        'BPix3+BPix4',
        'BPix1+FPix1_pos',
        'BPix1+FPix1_neg',
        'BPix2+FPix1_pos',
        'BPix2+FPix1_neg',
        'FPix1_pos+FPix2_pos',
        'FPix1_neg+FPix2_neg',
        'FPix1_pos+FPix3_pos',
        'FPix1_neg+FPix3_neg',
        'FPix2_pos+FPix3_pos',
        'FPix2_neg+FPix3_neg',
        'BPix4+TIB1',
        'BPix4+TIB2'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.jetCoreRegionalStepSeeds = cms.EDProducer("SeedCreatorFromRegionConsecutiveHitsEDProducer",
    MinOneOverPtError = cms.double(1),
    OriginTransverseErrorMultiplier = cms.double(1),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('none')
    ),
    SeedMomentumForBOFF = cms.double(5),
    TTRHBuilder = cms.string('WithTrackAngle'),
    forceKinematicWithRegionDirection = cms.bool(True),
    magneticField = cms.string('ParabolicMf'),
    mightGet = cms.optional.untracked.vstring,
    propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    seedingHitSets = cms.InputTag("jetCoreRegionalStepHitDoublets")
)


process.jetCoreRegionalStepSeedsBarrel = cms.EDProducer("DeepCoreSeedGenerator",
    centralMIPCharge = cms.double(2),
    chargeFractionMin = cms.double(18000),
    cores = cms.InputTag("jetsForCoreTrackingBarrel"),
    deltaR = cms.double(0.25),
    inputTensorName = cms.vstring(
        'input_1',
        'input_2',
        'input_3'
    ),
    mightGet = cms.optional.untracked.vstring,
    outputTensorName = cms.vstring(
        'output_node0',
        'output_node1'
    ),
    pixelCPE = cms.string('PixelCPEGeneric'),
    pixelClusters = cms.InputTag("siPixelClustersPreSplitting"),
    probThr = cms.double(0.85),
    ptMin = cms.double(100),
    vertices = cms.InputTag("firstStepPrimaryVertices"),
    weightFile = cms.FileInPath('RecoTracker/TkSeedGenerator/data/DeepCore/DeepCoreSeedGenerator_TrainedModel_barrel_2017.pb')
)


process.jetCoreRegionalStepSeedsEndcap = cms.EDProducer("SeedCreatorFromRegionConsecutiveHitsEDProducer",
    MinOneOverPtError = cms.double(1),
    OriginTransverseErrorMultiplier = cms.double(1),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('none')
    ),
    SeedMomentumForBOFF = cms.double(5),
    TTRHBuilder = cms.string('WithTrackAngle'),
    forceKinematicWithRegionDirection = cms.bool(True),
    magneticField = cms.string('ParabolicMf'),
    mightGet = cms.optional.untracked.vstring,
    propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    seedingHitSets = cms.InputTag("jetCoreRegionalStepEndcapHitDoublets")
)


process.jetCoreRegionalStepTrackCandidates = cms.EDProducer("CkfTrackCandidateMaker",
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    RedundantSeedCleaner = cms.string('CachingSeedCleanerBySharedInput'),
    TrajectoryBuilderPSet = cms.PSet(
        refToPSet_ = cms.string('jetCoreRegionalStepTrajectoryBuilder')
    ),
    TrajectoryCleaner = cms.string('TrajectoryCleanerBySharedHits'),
    TransientInitialStateEstimatorParameters = cms.PSet(
        numberMeasurementsForFit = cms.int32(4),
        propagatorAlongTISE = cms.string('PropagatorWithMaterial'),
        propagatorOppositeTISE = cms.string('PropagatorWithMaterialOpposite')
    ),
    cleanTrajectoryAfterInOut = cms.bool(True),
    clustersToSkip = cms.InputTag(""),
    doSeedingRegionRebuilding = cms.bool(True),
    maxNSeeds = cms.uint32(500000),
    maxSeedsBeforeCleaning = cms.uint32(10000),
    mightGet = cms.optional.untracked.vstring,
    numHitsForSeedCleaner = cms.int32(4),
    onlyPixelHitsForSeedCleaner = cms.bool(False),
    phase2clustersToSkip = cms.InputTag(""),
    reverseTrajectories = cms.bool(False),
    src = cms.InputTag("jetCoreRegionalStepSeeds"),
    useHitsSplitting = cms.bool(True)
)


process.jetCoreRegionalStepTrackingRegions = cms.EDProducer("TauRegionalPixelSeedTrackingRegionEDProducer",
    RegionPSet = cms.PSet(
        JetSrc = cms.InputTag("jetsForCoreTracking"),
        deltaEtaRegion = cms.double(0.2),
        deltaPhiRegion = cms.double(0.2),
        howToUseMeasurementTracker = cms.string('Never'),
        measurementTrackerName = cms.InputTag("MeasurementTrackerEvent"),
        originHalfLength = cms.double(0.2),
        originRadius = cms.double(0.2),
        ptMin = cms.double(10),
        searchOpt = cms.bool(False),
        vertexSrc = cms.InputTag("firstStepGoodPrimaryVertices")
    ),
    mightGet = cms.optional.untracked.vstring
)


process.jetCoreRegionalStepTracks = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('jetCoreRegionalStep'),
    Fitter = cms.string('FlexibleKFFittingSmoother'),
    GeometricInnerState = cms.bool(False),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    Propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    TTRHBuilder = cms.string('WithAngleAndTemplateWithoutProbQ'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("jetCoreRegionalStepTrackCandidates"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(True)
)


process.lowPtQuadStep = cms.EDProducer("TrackTfClassifier",
    beamspot = cms.InputTag("offlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mightGet = cms.optional.untracked.vstring,
    mva = cms.PSet(
        batchSize = cms.int32(16),
        tfDnnLabel = cms.string('trackSelectionTf_CKF')
    ),
    qualityCuts = cms.vdouble(-0.33, 0.13, 0.35),
    src = cms.InputTag("lowPtQuadStepTracks"),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.lowPtQuadStepClusters = cms.EDProducer("TrackClusterRemover",
    TrackQuality = cms.string('highPurity'),
    maxChi2 = cms.double(9.0),
    mightGet = cms.optional.untracked.vstring,
    minNumberOfLayersWithMeasBeforeFiltering = cms.int32(0),
    oldClusterRemovalInfo = cms.InputTag(""),
    overrideTrkQuals = cms.InputTag(""),
    pixelClusters = cms.InputTag("siPixelClusters"),
    stripClusters = cms.InputTag("siStripClusters"),
    trackClassifier = cms.InputTag("initialStep","QualityMasks"),
    trajectories = cms.InputTag("initialStepTracks")
)


process.lowPtQuadStepHitDoublets = cms.EDProducer("HitPairEDProducer",
    clusterCheck = cms.InputTag("trackerClusterCheck"),
    layerPairs = cms.vuint32(0, 1, 2),
    maxElement = cms.uint32(50000000),
    maxElementTotal = cms.uint32(50000000),
    mightGet = cms.optional.untracked.vstring,
    produceIntermediateHitDoublets = cms.bool(True),
    produceSeedingHitSets = cms.bool(False),
    seedingLayers = cms.InputTag("lowPtQuadStepSeedLayers"),
    trackingRegions = cms.InputTag("lowPtQuadStepTrackingRegions"),
    trackingRegionsSeedingLayers = cms.InputTag("")
)


process.lowPtQuadStepHitQuadruplets = cms.EDProducer("CAHitQuadrupletEDProducer",
    CAHardPtCut = cms.double(0),
    CAOnlyOneLastHitPerLayerFilter = cms.optional.bool,
    CAPhiCut = cms.double(0.3),
    CAPhiCut_byTriplets = cms.VPSet(cms.PSet(
        cut = cms.double(-1),
        seedingLayers = cms.string('')
    )),
    CAThetaCut = cms.double(0.0017),
    CAThetaCut_byTriplets = cms.VPSet(cms.PSet(
        cut = cms.double(-1),
        seedingLayers = cms.string('')
    )),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('LowPtClusterShapeSeedComparitor'),
        clusterShapeCacheSrc = cms.InputTag("siPixelClusterShapeCache"),
        clusterShapeHitFilter = cms.string('ClusterShapeHitFilter')
    ),
    doublets = cms.InputTag("lowPtQuadStepHitDoublets"),
    extraHitRPhitolerance = cms.double(0.032),
    fitFastCircle = cms.bool(True),
    fitFastCircleChi2Cut = cms.bool(True),
    maxChi2 = cms.PSet(
        enabled = cms.bool(True),
        pt1 = cms.double(0.7),
        pt2 = cms.double(2),
        value1 = cms.double(1000),
        value2 = cms.double(150)
    ),
    mightGet = cms.optional.untracked.vstring,
    useBendingCorrection = cms.bool(True)
)


process.lowPtQuadStepMasks = cms.EDProducer("FastTrackerRecHitMaskProducer",
    TrackQuality = cms.string('highPurity'),
    minNumberOfLayersWithMeasBeforeFiltering = cms.int32(0),
    minNumberOfLayersWithMeasBeforeFiltering_ = cms.int32(0),
    oldHitRemovalInfo = cms.InputTag(""),
    recHits = cms.InputTag("fastTrackerRecHits"),
    trackClassifier = cms.InputTag("initialStep","QualityMasks"),
    trackQuality = cms.string('hightPurity'),
    trajectories = cms.InputTag("initialStepTracks")
)


process.lowPtQuadStepSeedLayers = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(
        HitProducer = cms.string('siPixelRecHits'),
        TTRHBuilder = cms.string('WithTrackAngle'),
        skipClusters = cms.InputTag("lowPtQuadStepClusters")
    ),
    FPix = cms.PSet(
        HitProducer = cms.string('siPixelRecHits'),
        TTRHBuilder = cms.string('WithTrackAngle'),
        skipClusters = cms.InputTag("lowPtQuadStepClusters")
    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(

    ),
    TIB = cms.PSet(

    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(

    ),
    layerList = cms.vstring(
        'BPix1+BPix2+BPix3+BPix4',
        'BPix1+BPix2+BPix3+FPix1_pos',
        'BPix1+BPix2+BPix3+FPix1_neg',
        'BPix1+BPix2+FPix1_pos+FPix2_pos',
        'BPix1+BPix2+FPix1_neg+FPix2_neg',
        'BPix1+FPix1_pos+FPix2_pos+FPix3_pos',
        'BPix1+FPix1_neg+FPix2_neg+FPix3_neg'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.lowPtQuadStepSeeds = cms.EDProducer("SeedCreatorFromRegionConsecutiveHitsEDProducer",
    MinOneOverPtError = cms.double(1),
    OriginTransverseErrorMultiplier = cms.double(1),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('none')
    ),
    SeedMomentumForBOFF = cms.double(5),
    TTRHBuilder = cms.string('WithTrackAngle'),
    forceKinematicWithRegionDirection = cms.bool(False),
    magneticField = cms.string('ParabolicMf'),
    mightGet = cms.optional.untracked.vstring,
    propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    seedingHitSets = cms.InputTag("lowPtQuadStepHitQuadruplets")
)


process.lowPtQuadStepSelector = cms.EDProducer("MultiTrackSelector",
    beamspot = cms.InputTag("offlineBeamSpot"),
    src = cms.InputTag("lowPtQuadStepTracks"),
    trackSelectors = cms.VPSet(
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(True),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(2.0),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(0.8, 4.0),
            d0_par2 = cms.vdouble(0.6, 4.0),
            dz_par1 = cms.vdouble(0.7, 4.0),
            dz_par2 = cms.vdouble(0.6, 4.0),
            keepAllTracks = cms.bool(False),
            maxNumberLostLayers = cms.uint32(2),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(3),
            minNumberLayers = cms.uint32(3),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(0),
            nSigmaZ = cms.double(4.0),
            name = cms.string('lowPtQuadStepLoose'),
            preFilterName = cms.string(''),
            qualityBit = cms.string('loose'),
            res_par = cms.vdouble(0.003, 0.002),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        ),
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(True),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(1.4),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(0.7, 4.0),
            d0_par2 = cms.vdouble(0.5, 4.0),
            dz_par1 = cms.vdouble(0.6, 4.0),
            dz_par2 = cms.vdouble(0.5, 4.0),
            keepAllTracks = cms.bool(True),
            maxNumberLostLayers = cms.uint32(2),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(3),
            minNumberLayers = cms.uint32(3),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(0),
            nSigmaZ = cms.double(4.0),
            name = cms.string('lowPtQuadStepTight'),
            preFilterName = cms.string('lowPtQuadStepLoose'),
            qualityBit = cms.string('tight'),
            res_par = cms.vdouble(0.003, 0.002),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        ),
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(True),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(1.2),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(0.5, 4.0),
            d0_par2 = cms.vdouble(0.45, 4.0),
            dz_par1 = cms.vdouble(0.5, 4.0),
            dz_par2 = cms.vdouble(0.45, 4.0),
            keepAllTracks = cms.bool(True),
            maxNumberLostLayers = cms.uint32(2),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(4.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(3),
            minNumberLayers = cms.uint32(3),
            min_eta = cms.double(-4.0),
            min_nhits = cms.uint32(0),
            nSigmaZ = cms.double(4.0),
            name = cms.string('lowPtQuadStep'),
            preFilterName = cms.string('lowPtQuadStepTight'),
            qualityBit = cms.string('highPurity'),
            res_par = cms.vdouble(0.003, 0.001),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        )
    ),
    useVertices = cms.bool(True),
    useVtxError = cms.bool(False),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.lowPtQuadStepTrackCandidates = cms.EDProducer("CkfTrackCandidateMaker",
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    RedundantSeedCleaner = cms.string('CachingSeedCleanerBySharedInput'),
    TrajectoryBuilderPSet = cms.PSet(
        refToPSet_ = cms.string('lowPtQuadStepTrajectoryBuilder')
    ),
    TrajectoryCleaner = cms.string('lowPtQuadStepTrajectoryCleanerBySharedHits'),
    TransientInitialStateEstimatorParameters = cms.PSet(
        numberMeasurementsForFit = cms.int32(4),
        propagatorAlongTISE = cms.string('PropagatorWithMaterial'),
        propagatorOppositeTISE = cms.string('PropagatorWithMaterialOpposite')
    ),
    cleanTrajectoryAfterInOut = cms.bool(True),
    clustersToSkip = cms.InputTag("lowPtQuadStepClusters"),
    doSeedingRegionRebuilding = cms.bool(True),
    maxNSeeds = cms.uint32(500000),
    maxSeedsBeforeCleaning = cms.uint32(5000),
    mightGet = cms.optional.untracked.vstring,
    numHitsForSeedCleaner = cms.int32(50),
    onlyPixelHitsForSeedCleaner = cms.bool(True),
    phase2clustersToSkip = cms.InputTag(""),
    reverseTrajectories = cms.bool(False),
    src = cms.InputTag("lowPtQuadStepSeeds"),
    useHitsSplitting = cms.bool(True)
)


process.lowPtQuadStepTrackCandidatesMkFit = cms.EDProducer("MkFitProducer",
    backwardFitInCMSSW = cms.bool(False),
    buildingRoutine = cms.string('cloneEngine'),
    clustersToSkip = cms.InputTag("lowPtQuadStepClusters"),
    config = cms.ESInputTag("","lowPtQuadStepTrackCandidatesMkFitConfig"),
    eventOfHits = cms.InputTag("mkFitEventOfHits"),
    limitConcurrency = cms.untracked.bool(False),
    mightGet = cms.optional.untracked.vstring,
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutLoose')
    ),
    mkFitSilent = cms.untracked.bool(True),
    pixelHits = cms.InputTag("mkFitSiPixelHits"),
    removeDuplicates = cms.bool(True),
    seedCleaning = cms.bool(True),
    seeds = cms.InputTag("lowPtQuadStepTrackCandidatesMkFitSeeds"),
    stripHits = cms.InputTag("mkFitSiStripHits")
)


process.lowPtQuadStepTrackCandidatesMkFitSeeds = cms.EDProducer("MkFitSeedConverter",
    maxNSeeds = cms.uint32(500000),
    mightGet = cms.optional.untracked.vstring,
    seeds = cms.InputTag("lowPtQuadStepSeeds"),
    ttrhBuilder = cms.ESInputTag("","WithTrackAngle")
)


process.lowPtQuadStepTrackingRegions = cms.EDProducer("GlobalTrackingRegionFromBeamSpotEDProducer",
    RegionPSet = cms.PSet(
        beamSpot = cms.InputTag("offlineBeamSpot"),
        nSigmaZ = cms.double(4),
        originHalfLength = cms.double(0),
        originRadius = cms.double(0.02),
        precise = cms.bool(True),
        ptMin = cms.double(0.15),
        useMultipleScattering = cms.bool(False)
    ),
    mightGet = cms.optional.untracked.vstring
)


process.lowPtQuadStepTracks = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('lowPtQuadStep'),
    Fitter = cms.string('FlexibleKFFittingSmoother'),
    GeometricInnerState = cms.bool(False),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    Propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    TTRHBuilder = cms.string('WithAngleAndTemplateWithoutProbQ'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("lowPtQuadStepTrackCandidates"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(True)
)


process.lowPtTripletStep = cms.EDProducer("TrackTfClassifier",
    beamspot = cms.InputTag("offlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mightGet = cms.optional.untracked.vstring,
    mva = cms.PSet(
        batchSize = cms.int32(16),
        tfDnnLabel = cms.string('trackSelectionTf_CKF')
    ),
    qualityCuts = cms.vdouble(-0.23, 0.15, 0.41),
    src = cms.InputTag("lowPtTripletStepTracks"),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.lowPtTripletStepClusters = cms.EDProducer("TrackClusterRemover",
    TrackQuality = cms.string('highPurity'),
    maxChi2 = cms.double(9.0),
    mightGet = cms.optional.untracked.vstring,
    minNumberOfLayersWithMeasBeforeFiltering = cms.int32(0),
    oldClusterRemovalInfo = cms.InputTag("highPtTripletStepClusters"),
    overrideTrkQuals = cms.InputTag(""),
    pixelClusters = cms.InputTag("siPixelClusters"),
    stripClusters = cms.InputTag("siStripClusters"),
    trackClassifier = cms.InputTag("highPtTripletStep","QualityMasks"),
    trajectories = cms.InputTag("highPtTripletStepTracks")
)


process.lowPtTripletStepHitDoublets = cms.EDProducer("HitPairEDProducer",
    clusterCheck = cms.InputTag("trackerClusterCheck"),
    layerPairs = cms.vuint32(0, 1),
    maxElement = cms.uint32(50000000),
    maxElementTotal = cms.uint32(50000000),
    mightGet = cms.optional.untracked.vstring,
    produceIntermediateHitDoublets = cms.bool(True),
    produceSeedingHitSets = cms.bool(False),
    seedingLayers = cms.InputTag("lowPtTripletStepSeedLayers"),
    trackingRegions = cms.InputTag("lowPtTripletStepTrackingRegions"),
    trackingRegionsSeedingLayers = cms.InputTag("")
)


process.lowPtTripletStepHitTriplets = cms.EDProducer("CAHitTripletEDProducer",
    CAHardPtCut = cms.double(0),
    CAPhiCut = cms.double(0.05),
    CAPhiCut_byTriplets = cms.VPSet(cms.PSet(
        cut = cms.double(-1),
        seedingLayers = cms.string('')
    )),
    CAThetaCut = cms.double(0.002),
    CAThetaCut_byTriplets = cms.VPSet(cms.PSet(
        cut = cms.double(-1),
        seedingLayers = cms.string('')
    )),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('LowPtClusterShapeSeedComparitor'),
        clusterShapeCacheSrc = cms.InputTag("siPixelClusterShapeCache"),
        clusterShapeHitFilter = cms.string('ClusterShapeHitFilter')
    ),
    doublets = cms.InputTag("lowPtTripletStepHitDoublets"),
    extraHitRPhitolerance = cms.double(0.032),
    maxChi2 = cms.PSet(
        enabled = cms.bool(True),
        pt1 = cms.double(0.8),
        pt2 = cms.double(2),
        value1 = cms.double(70),
        value2 = cms.double(8)
    ),
    mightGet = cms.optional.untracked.vstring,
    useBendingCorrection = cms.bool(True)
)


process.lowPtTripletStepMasks = cms.EDProducer("FastTrackerRecHitMaskProducer",
    TrackQuality = cms.string('highPurity'),
    minNumberOfLayersWithMeasBeforeFiltering = cms.int32(0),
    minNumberOfLayersWithMeasBeforeFiltering_ = cms.int32(0),
    oldHitRemovalInfo = cms.InputTag("highPtTripletStepMasks"),
    recHits = cms.InputTag("fastTrackerRecHits"),
    trackClassifier = cms.InputTag("highPtTripletStep","QualityMasks"),
    trackQuality = cms.string('hightPurity'),
    trajectories = cms.InputTag("highPtTripletStepTracks")
)


process.lowPtTripletStepSeedLayers = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(
        HitProducer = cms.string('siPixelRecHits'),
        TTRHBuilder = cms.string('WithTrackAngle'),
        skipClusters = cms.InputTag("lowPtTripletStepClusters")
    ),
    FPix = cms.PSet(
        HitProducer = cms.string('siPixelRecHits'),
        TTRHBuilder = cms.string('WithTrackAngle'),
        skipClusters = cms.InputTag("lowPtTripletStepClusters")
    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(

    ),
    TIB = cms.PSet(

    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(

    ),
    layerList = cms.vstring(
        'BPix1+BPix2+BPix3',
        'BPix2+BPix3+BPix4',
        'BPix1+BPix3+BPix4',
        'BPix1+BPix2+BPix4',
        'BPix2+BPix3+FPix1_pos',
        'BPix2+BPix3+FPix1_neg',
        'BPix1+BPix2+FPix1_pos',
        'BPix1+BPix2+FPix1_neg',
        'BPix1+BPix3+FPix1_pos',
        'BPix1+BPix3+FPix1_neg',
        'BPix2+FPix1_pos+FPix2_pos',
        'BPix2+FPix1_neg+FPix2_neg',
        'BPix1+FPix1_pos+FPix2_pos',
        'BPix1+FPix1_neg+FPix2_neg',
        'BPix1+BPix2+FPix2_pos',
        'BPix1+BPix2+FPix2_neg',
        'FPix1_pos+FPix2_pos+FPix3_pos',
        'FPix1_neg+FPix2_neg+FPix3_neg',
        'BPix1+FPix2_pos+FPix3_pos',
        'BPix1+FPix2_neg+FPix3_neg',
        'BPix1+FPix1_pos+FPix3_pos',
        'BPix1+FPix1_neg+FPix3_neg'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.lowPtTripletStepSeeds = cms.EDProducer("SeedCreatorFromRegionConsecutiveHitsEDProducer",
    MinOneOverPtError = cms.double(1),
    OriginTransverseErrorMultiplier = cms.double(1),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('none')
    ),
    SeedMomentumForBOFF = cms.double(5),
    TTRHBuilder = cms.string('WithTrackAngle'),
    forceKinematicWithRegionDirection = cms.bool(False),
    magneticField = cms.string('ParabolicMf'),
    mightGet = cms.optional.untracked.vstring,
    propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    seedingHitSets = cms.InputTag("lowPtTripletStepHitTriplets")
)


process.lowPtTripletStepSelector = cms.EDProducer("MultiTrackSelector",
    GBRForestLabel = cms.string('MVASelectorIter1'),
    beamspot = cms.InputTag("offlineBeamSpot"),
    src = cms.InputTag("lowPtTripletStepTracks"),
    trackSelectors = cms.VPSet(
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(True),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(1.6),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(0.55, 4.0),
            d0_par2 = cms.vdouble(0.55, 4.0),
            dz_par1 = cms.vdouble(0.65, 4.0),
            dz_par2 = cms.vdouble(0.45, 4.0),
            keepAllTracks = cms.bool(False),
            maxNumberLostLayers = cms.uint32(999),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(0),
            minNumberLayers = cms.uint32(0),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(0),
            nSigmaZ = cms.double(4.0),
            name = cms.string('lowPtTripletStepLoose'),
            preFilterName = cms.string(''),
            qualityBit = cms.string('loose'),
            res_par = cms.vdouble(0.003, 0.01),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        ),
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(True),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(0.7),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(0.3, 4.0),
            d0_par2 = cms.vdouble(0.4, 4.0),
            dz_par1 = cms.vdouble(0.35, 4.0),
            dz_par2 = cms.vdouble(0.4, 4.0),
            keepAllTracks = cms.bool(True),
            maxNumberLostLayers = cms.uint32(2),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(3),
            minNumberLayers = cms.uint32(3),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(0),
            nSigmaZ = cms.double(4.0),
            name = cms.string('lowPtTripletStepTight'),
            preFilterName = cms.string('lowPtTripletStepLoose'),
            qualityBit = cms.string('tight'),
            res_par = cms.vdouble(0.003, 0.01),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        ),
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(True),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(0.7),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(0.3, 4.0),
            d0_par2 = cms.vdouble(0.4, 4.0),
            dz_par1 = cms.vdouble(0.35, 4.0),
            dz_par2 = cms.vdouble(0.4, 4.0),
            keepAllTracks = cms.bool(True),
            maxNumberLostLayers = cms.uint32(2),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(3),
            minNumberLayers = cms.uint32(3),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(0),
            nSigmaZ = cms.double(4.0),
            name = cms.string('QualityMasks'),
            preFilterName = cms.string('lowPtTripletStepTight'),
            qualityBit = cms.string('highPurity'),
            res_par = cms.vdouble(0.003, 0.001),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        )
    ),
    useAnyMVA = cms.bool(False),
    useVertices = cms.bool(True),
    useVtxError = cms.bool(False),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.lowPtTripletStepTrackCandidates = cms.EDProducer("CkfTrackCandidateMaker",
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    RedundantSeedCleaner = cms.string('CachingSeedCleanerBySharedInput'),
    TrajectoryBuilderPSet = cms.PSet(
        refToPSet_ = cms.string('lowPtTripletStepTrajectoryBuilder')
    ),
    TrajectoryCleaner = cms.string('lowPtTripletStepTrajectoryCleanerBySharedHits'),
    TransientInitialStateEstimatorParameters = cms.PSet(
        numberMeasurementsForFit = cms.int32(4),
        propagatorAlongTISE = cms.string('PropagatorWithMaterial'),
        propagatorOppositeTISE = cms.string('PropagatorWithMaterialOpposite')
    ),
    cleanTrajectoryAfterInOut = cms.bool(True),
    clustersToSkip = cms.InputTag("lowPtTripletStepClusters"),
    doSeedingRegionRebuilding = cms.bool(True),
    maxNSeeds = cms.uint32(500000),
    maxSeedsBeforeCleaning = cms.uint32(5000),
    mightGet = cms.optional.untracked.vstring,
    numHitsForSeedCleaner = cms.int32(50),
    onlyPixelHitsForSeedCleaner = cms.bool(True),
    phase2clustersToSkip = cms.InputTag(""),
    reverseTrajectories = cms.bool(False),
    src = cms.InputTag("lowPtTripletStepSeeds"),
    useHitsSplitting = cms.bool(True)
)


process.lowPtTripletStepTrackCandidatesMkFit = cms.EDProducer("MkFitProducer",
    backwardFitInCMSSW = cms.bool(False),
    buildingRoutine = cms.string('cloneEngine'),
    clustersToSkip = cms.InputTag("lowPtTripletStepClusters"),
    config = cms.ESInputTag("","lowPtTripletStepTrackCandidatesMkFitConfig"),
    eventOfHits = cms.InputTag("mkFitEventOfHits"),
    limitConcurrency = cms.untracked.bool(False),
    mightGet = cms.optional.untracked.vstring,
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutLoose')
    ),
    mkFitSilent = cms.untracked.bool(True),
    pixelHits = cms.InputTag("mkFitSiPixelHits"),
    removeDuplicates = cms.bool(True),
    seedCleaning = cms.bool(True),
    seeds = cms.InputTag("lowPtTripletStepTrackCandidatesMkFitSeeds"),
    stripHits = cms.InputTag("mkFitSiStripHits")
)


process.lowPtTripletStepTrackCandidatesMkFitSeeds = cms.EDProducer("MkFitSeedConverter",
    maxNSeeds = cms.uint32(500000),
    mightGet = cms.optional.untracked.vstring,
    seeds = cms.InputTag("lowPtTripletStepSeeds"),
    ttrhBuilder = cms.ESInputTag("","WithTrackAngle")
)


process.lowPtTripletStepTrackingRegions = cms.EDProducer("GlobalTrackingRegionFromBeamSpotEDProducer",
    RegionPSet = cms.PSet(
        beamSpot = cms.InputTag("offlineBeamSpot"),
        nSigmaZ = cms.double(4),
        originHalfLength = cms.double(0),
        originRadius = cms.double(0.02),
        precise = cms.bool(True),
        ptMin = cms.double(0.2),
        useMultipleScattering = cms.bool(False)
    ),
    mightGet = cms.optional.untracked.vstring
)


process.lowPtTripletStepTracks = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('lowPtTripletStep'),
    Fitter = cms.string('FlexibleKFFittingSmoother'),
    GeometricInnerState = cms.bool(False),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    Propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    TTRHBuilder = cms.string('WithAngleAndTemplateWithoutProbQ'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("lowPtTripletStepTrackCandidates"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(True)
)


process.mergedDuplicateTracks = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('undefAlgorithm'),
    Fitter = cms.string('RKFittingSmoother'),
    GeometricInnerState = cms.bool(False),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    Propagator = cms.string('RungeKuttaTrackerPropagator'),
    SimpleMagneticField = cms.string(''),
    TTRHBuilder = cms.string('WithAngleAndTemplate'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("duplicateTrackCandidates","candidates"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(False)
)


process.mixedTripletStep = cms.EDProducer("TrackTfClassifier",
    beamspot = cms.InputTag("offlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mightGet = cms.optional.untracked.vstring,
    mva = cms.PSet(
        batchSize = cms.int32(16),
        tfDnnLabel = cms.string('trackSelectionTf')
    ),
    qualityCuts = cms.vdouble(-0.86, -0.68, -0.43),
    src = cms.InputTag("mixedTripletStepTracks"),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.mixedTripletStepClassifier1 = cms.EDProducer("TrackMVAClassifierDetached",
    beamspot = cms.InputTag("offlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mightGet = cms.optional.untracked.vstring,
    mva = cms.PSet(
        GBRForestFileName = cms.string(''),
        GBRForestLabel = cms.string('MVASelectorIter4_13TeV')
    ),
    qualityCuts = cms.vdouble(-0.5, 0.0, 0.5),
    src = cms.InputTag("mixedTripletStepTracks"),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.mixedTripletStepClassifier2 = cms.EDProducer("TrackMVAClassifierPrompt",
    beamspot = cms.InputTag("offlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mightGet = cms.optional.untracked.vstring,
    mva = cms.PSet(
        GBRForestFileName = cms.string(''),
        GBRForestLabel = cms.string('MVASelectorIter0_13TeV')
    ),
    qualityCuts = cms.vdouble(-0.2, -0.2, -0.2),
    src = cms.InputTag("mixedTripletStepTracks"),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.mixedTripletStepClusters = cms.EDProducer("TrackClusterRemover",
    TrackQuality = cms.string('highPurity'),
    maxChi2 = cms.double(9.0),
    mightGet = cms.optional.untracked.vstring,
    minNumberOfLayersWithMeasBeforeFiltering = cms.int32(0),
    oldClusterRemovalInfo = cms.InputTag("chargeCut2069Clusters"),
    overrideTrkQuals = cms.InputTag(""),
    pixelClusters = cms.InputTag("siPixelClusters"),
    stripClusters = cms.InputTag("siStripClusters"),
    trackClassifier = cms.InputTag("pixelPairStep","QualityMasks"),
    trajectories = cms.InputTag("pixelPairStepTracks")
)


process.mixedTripletStepHitDoubletsA = cms.EDProducer("HitPairEDProducer",
    clusterCheck = cms.InputTag("trackerClusterCheck"),
    layerPairs = cms.vuint32(0),
    maxElement = cms.uint32(50000000),
    maxElementTotal = cms.uint32(50000000),
    mightGet = cms.optional.untracked.vstring,
    produceIntermediateHitDoublets = cms.bool(True),
    produceSeedingHitSets = cms.bool(False),
    seedingLayers = cms.InputTag("mixedTripletStepSeedLayersA"),
    trackingRegions = cms.InputTag("mixedTripletStepTrackingRegionsA"),
    trackingRegionsSeedingLayers = cms.InputTag("")
)


process.mixedTripletStepHitDoubletsB = cms.EDProducer("HitPairEDProducer",
    clusterCheck = cms.InputTag("trackerClusterCheck"),
    layerPairs = cms.vuint32(0),
    maxElement = cms.uint32(50000000),
    maxElementTotal = cms.uint32(50000000),
    mightGet = cms.optional.untracked.vstring,
    produceIntermediateHitDoublets = cms.bool(True),
    produceSeedingHitSets = cms.bool(False),
    seedingLayers = cms.InputTag("mixedTripletStepSeedLayersB"),
    trackingRegions = cms.InputTag("mixedTripletStepTrackingRegionsB"),
    trackingRegionsSeedingLayers = cms.InputTag("")
)


process.mixedTripletStepHitTripletsA = cms.EDProducer("PixelTripletLargeTipEDProducer",
    doublets = cms.InputTag("mixedTripletStepHitDoubletsA"),
    extraHitRPhitolerance = cms.double(0),
    extraHitRZtolerance = cms.double(0),
    maxElement = cms.uint32(1000000),
    mightGet = cms.optional.untracked.vstring,
    phiPreFiltering = cms.double(0.3),
    produceIntermediateHitTriplets = cms.bool(False),
    produceSeedingHitSets = cms.bool(True),
    useBending = cms.bool(True),
    useFixedPreFiltering = cms.bool(False),
    useMultScattering = cms.bool(True)
)


process.mixedTripletStepHitTripletsB = cms.EDProducer("PixelTripletLargeTipEDProducer",
    doublets = cms.InputTag("mixedTripletStepHitDoubletsB"),
    extraHitRPhitolerance = cms.double(0),
    extraHitRZtolerance = cms.double(0),
    maxElement = cms.uint32(1000000),
    mightGet = cms.optional.untracked.vstring,
    phiPreFiltering = cms.double(0.3),
    produceIntermediateHitTriplets = cms.bool(False),
    produceSeedingHitSets = cms.bool(True),
    useBending = cms.bool(True),
    useFixedPreFiltering = cms.bool(False),
    useMultScattering = cms.bool(True)
)


process.mixedTripletStepMasks = cms.EDProducer("FastTrackerRecHitMaskProducer",
    TrackQuality = cms.string('highPurity'),
    minNumberOfLayersWithMeasBeforeFiltering = cms.int32(0),
    minNumberOfLayersWithMeasBeforeFiltering_ = cms.int32(0),
    oldHitRemovalInfo = cms.InputTag("pixelPairStepMasks"),
    recHits = cms.InputTag("fastTrackerRecHits"),
    trackClassifier = cms.InputTag("pixelPairStep","QualityMasks"),
    trackQuality = cms.string('hightPurity'),
    trajectories = cms.InputTag("pixelPairStepTracks")
)


process.mixedTripletStepSeedClusterMask = cms.EDProducer("SeedClusterRemover",
    Common = cms.PSet(
        maxChi2 = cms.double(9.0)
    ),
    oldClusterRemovalInfo = cms.InputTag("detachedTripletStepSeedClusterMask"),
    pixelClusters = cms.InputTag("siPixelClusters"),
    stripClusters = cms.InputTag("siStripClusters"),
    trajectories = cms.InputTag("mixedTripletStepSeeds")
)


process.mixedTripletStepSeedLayersA = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(
        HitProducer = cms.string('siPixelRecHits'),
        TTRHBuilder = cms.string('WithTrackAngle'),
        skipClusters = cms.InputTag("mixedTripletStepClusters")
    ),
    FPix = cms.PSet(
        HitProducer = cms.string('siPixelRecHits'),
        TTRHBuilder = cms.string('WithTrackAngle'),
        skipClusters = cms.InputTag("mixedTripletStepClusters")
    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTight')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        maxRing = cms.int32(1),
        minRing = cms.int32(1),
        skipClusters = cms.InputTag("mixedTripletStepClusters"),
        useRingSlector = cms.bool(True)
    ),
    TIB = cms.PSet(

    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(

    ),
    layerList = cms.vstring(
        'BPix2+FPix1_pos+FPix2_pos',
        'BPix2+FPix1_neg+FPix2_neg'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.mixedTripletStepSeedLayersB = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(
        HitProducer = cms.string('siPixelRecHits'),
        TTRHBuilder = cms.string('WithTrackAngle'),
        skipClusters = cms.InputTag("mixedTripletStepClusters")
    ),
    FPix = cms.PSet(

    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(

    ),
    TIB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTight')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        skipClusters = cms.InputTag("mixedTripletStepClusters")
    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(

    ),
    layerList = cms.vstring('BPix3+BPix4+TIB1'),
    mightGet = cms.optional.untracked.vstring
)


process.mixedTripletStepSeeds = cms.EDProducer("SeedCombiner",
    seedCollections = cms.VInputTag("mixedTripletStepSeedsA", "mixedTripletStepSeedsB")
)


process.mixedTripletStepSeedsA = cms.EDProducer("SeedCreatorFromRegionConsecutiveHitsTripletOnlyEDProducer",
    MinOneOverPtError = cms.double(1),
    OriginTransverseErrorMultiplier = cms.double(1),
    SeedComparitorPSet = cms.PSet(
        ClusterShapeCacheSrc = cms.InputTag("siPixelClusterShapeCache"),
        ClusterShapeHitFilterName = cms.string('mixedTripletStepClusterShapeHitFilter'),
        ComponentName = cms.string('PixelClusterShapeSeedComparitor'),
        FilterAtHelixStage = cms.bool(False),
        FilterPixelHits = cms.bool(True),
        FilterStripHits = cms.bool(True)
    ),
    SeedMomentumForBOFF = cms.double(5),
    TTRHBuilder = cms.string('WithTrackAngle'),
    forceKinematicWithRegionDirection = cms.bool(False),
    magneticField = cms.string('ParabolicMf'),
    mightGet = cms.optional.untracked.vstring,
    propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    seedingHitSets = cms.InputTag("mixedTripletStepHitTripletsA")
)


process.mixedTripletStepSeedsB = cms.EDProducer("SeedCreatorFromRegionConsecutiveHitsTripletOnlyEDProducer",
    MinOneOverPtError = cms.double(1),
    OriginTransverseErrorMultiplier = cms.double(1),
    SeedComparitorPSet = cms.PSet(
        ClusterShapeCacheSrc = cms.InputTag("siPixelClusterShapeCache"),
        ClusterShapeHitFilterName = cms.string('mixedTripletStepClusterShapeHitFilter'),
        ComponentName = cms.string('PixelClusterShapeSeedComparitor'),
        FilterAtHelixStage = cms.bool(False),
        FilterPixelHits = cms.bool(True),
        FilterStripHits = cms.bool(True)
    ),
    SeedMomentumForBOFF = cms.double(5),
    TTRHBuilder = cms.string('WithTrackAngle'),
    forceKinematicWithRegionDirection = cms.bool(False),
    magneticField = cms.string('ParabolicMf'),
    mightGet = cms.optional.untracked.vstring,
    propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    seedingHitSets = cms.InputTag("mixedTripletStepHitTripletsB")
)


process.mixedTripletStepSelector = cms.EDProducer("MultiTrackSelector",
    GBRForestLabel = cms.string('MVASelectorIter4'),
    beamspot = cms.InputTag("offlineBeamSpot"),
    src = cms.InputTag("mixedTripletStepTracks"),
    trackSelectors = cms.VPSet(
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(True),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(1.2),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(1.2, 3.0),
            d0_par2 = cms.vdouble(1.3, 3.0),
            dz_par1 = cms.vdouble(1.2, 3.0),
            dz_par2 = cms.vdouble(1.3, 3.0),
            keepAllTracks = cms.bool(False),
            maxNumberLostLayers = cms.uint32(1),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(2),
            minNumberLayers = cms.uint32(3),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(0),
            nSigmaZ = cms.double(4.0),
            name = cms.string('mixedTripletStepVtxLoose'),
            preFilterName = cms.string(''),
            qualityBit = cms.string('loose'),
            res_par = cms.vdouble(0.003, 0.001),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        ),
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(True),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(0.6),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(1.2, 4.0),
            d0_par2 = cms.vdouble(1.2, 4.0),
            dz_par1 = cms.vdouble(1.2, 4.0),
            dz_par2 = cms.vdouble(1.2, 4.0),
            keepAllTracks = cms.bool(False),
            maxNumberLostLayers = cms.uint32(1),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(3),
            minNumberLayers = cms.uint32(4),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(0),
            nSigmaZ = cms.double(4.0),
            name = cms.string('mixedTripletStepTrkLoose'),
            preFilterName = cms.string(''),
            qualityBit = cms.string('loose'),
            res_par = cms.vdouble(0.003, 0.001),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        ),
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(True),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(0.6),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(1.1, 3.0),
            d0_par2 = cms.vdouble(1.2, 3.0),
            dz_par1 = cms.vdouble(1.1, 3.0),
            dz_par2 = cms.vdouble(1.2, 3.0),
            keepAllTracks = cms.bool(True),
            maxNumberLostLayers = cms.uint32(1),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(3),
            minNumberLayers = cms.uint32(3),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(0),
            nSigmaZ = cms.double(4.0),
            name = cms.string('mixedTripletStepVtxTight'),
            preFilterName = cms.string('mixedTripletStepVtxLoose'),
            qualityBit = cms.string('tight'),
            res_par = cms.vdouble(0.003, 0.001),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        ),
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(True),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(0.4),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(1.1, 4.0),
            d0_par2 = cms.vdouble(1.1, 4.0),
            dz_par1 = cms.vdouble(1.1, 4.0),
            dz_par2 = cms.vdouble(1.1, 4.0),
            keepAllTracks = cms.bool(True),
            maxNumberLostLayers = cms.uint32(1),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(4),
            minNumberLayers = cms.uint32(5),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(0),
            nSigmaZ = cms.double(4.0),
            name = cms.string('mixedTripletStepTrkTight'),
            preFilterName = cms.string('mixedTripletStepTrkLoose'),
            qualityBit = cms.string('tight'),
            res_par = cms.vdouble(0.003, 0.001),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        ),
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(True),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(0.4),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(1.1, 3.0),
            d0_par2 = cms.vdouble(1.2, 3.0),
            dz_par1 = cms.vdouble(1.1, 3.0),
            dz_par2 = cms.vdouble(1.2, 3.0),
            keepAllTracks = cms.bool(True),
            maxNumberLostLayers = cms.uint32(1),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(3),
            minNumberLayers = cms.uint32(3),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(0),
            nSigmaZ = cms.double(4.0),
            name = cms.string('mixedTripletStepVtx'),
            preFilterName = cms.string('mixedTripletStepVtxTight'),
            qualityBit = cms.string('highPurity'),
            res_par = cms.vdouble(0.003, 0.001),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        ),
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(True),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(0.3),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(0.9, 4.0),
            d0_par2 = cms.vdouble(0.9, 4.0),
            dz_par1 = cms.vdouble(0.9, 4.0),
            dz_par2 = cms.vdouble(0.9, 4.0),
            keepAllTracks = cms.bool(True),
            maxNumberLostLayers = cms.uint32(0),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(4),
            minNumberLayers = cms.uint32(5),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(0),
            nSigmaZ = cms.double(4.0),
            name = cms.string('mixedTripletStepTrk'),
            preFilterName = cms.string('mixedTripletStepTrkTight'),
            qualityBit = cms.string('highPurity'),
            res_par = cms.vdouble(0.003, 0.001),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        )
    ),
    useAnyMVA = cms.bool(False),
    useVertices = cms.bool(True),
    useVtxError = cms.bool(False),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.mixedTripletStepTrackCandidates = cms.EDProducer("CkfTrackCandidateMaker",
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    RedundantSeedCleaner = cms.string('CachingSeedCleanerBySharedInput'),
    TrajectoryBuilderPSet = cms.PSet(
        refToPSet_ = cms.string('mixedTripletStepTrajectoryBuilder')
    ),
    TrajectoryCleaner = cms.string('mixedTripletStepTrajectoryCleanerBySharedHits'),
    TransientInitialStateEstimatorParameters = cms.PSet(
        numberMeasurementsForFit = cms.int32(4),
        propagatorAlongTISE = cms.string('PropagatorWithMaterial'),
        propagatorOppositeTISE = cms.string('PropagatorWithMaterialOpposite')
    ),
    cleanTrajectoryAfterInOut = cms.bool(True),
    clustersToSkip = cms.InputTag("mixedTripletStepClusters"),
    doSeedingRegionRebuilding = cms.bool(True),
    maxNSeeds = cms.uint32(500000),
    maxSeedsBeforeCleaning = cms.uint32(5000),
    mightGet = cms.optional.untracked.vstring,
    numHitsForSeedCleaner = cms.int32(50),
    onlyPixelHitsForSeedCleaner = cms.bool(False),
    phase2clustersToSkip = cms.InputTag(""),
    reverseTrajectories = cms.bool(False),
    src = cms.InputTag("mixedTripletStepSeeds"),
    useHitsSplitting = cms.bool(True)
)


process.mixedTripletStepTrackCandidatesMkFit = cms.EDProducer("MkFitProducer",
    backwardFitInCMSSW = cms.bool(False),
    buildingRoutine = cms.string('cloneEngine'),
    clustersToSkip = cms.InputTag("mixedTripletStepClusters"),
    config = cms.ESInputTag("","mixedTripletStepTrackCandidatesMkFitConfig"),
    eventOfHits = cms.InputTag("mkFitEventOfHits"),
    limitConcurrency = cms.untracked.bool(False),
    mightGet = cms.optional.untracked.vstring,
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutLoose')
    ),
    mkFitSilent = cms.untracked.bool(True),
    pixelHits = cms.InputTag("mkFitSiPixelHits"),
    removeDuplicates = cms.bool(True),
    seedCleaning = cms.bool(True),
    seeds = cms.InputTag("mixedTripletStepTrackCandidatesMkFitSeeds"),
    stripHits = cms.InputTag("mkFitSiStripHits")
)


process.mixedTripletStepTrackCandidatesMkFitSeeds = cms.EDProducer("MkFitSeedConverter",
    maxNSeeds = cms.uint32(500000),
    mightGet = cms.optional.untracked.vstring,
    seeds = cms.InputTag("mixedTripletStepSeeds"),
    ttrhBuilder = cms.ESInputTag("","WithTrackAngle")
)


process.mixedTripletStepTrackingRegionsA = cms.EDProducer("GlobalTrackingRegionFromBeamSpotEDProducer",
    RegionPSet = cms.PSet(
        beamSpot = cms.InputTag("offlineBeamSpot"),
        nSigmaZ = cms.double(0),
        originHalfLength = cms.double(15.0),
        originRadius = cms.double(1.5),
        precise = cms.bool(True),
        ptMin = cms.double(0.4),
        useMultipleScattering = cms.bool(False)
    ),
    mightGet = cms.optional.untracked.vstring
)


process.mixedTripletStepTrackingRegionsB = cms.EDProducer("GlobalTrackingRegionFromBeamSpotEDProducer",
    RegionPSet = cms.PSet(
        beamSpot = cms.InputTag("offlineBeamSpot"),
        nSigmaZ = cms.double(0),
        originHalfLength = cms.double(10.0),
        originRadius = cms.double(1.5),
        precise = cms.bool(True),
        ptMin = cms.double(0.6),
        useMultipleScattering = cms.bool(False)
    ),
    mightGet = cms.optional.untracked.vstring
)


process.mixedTripletStepTracks = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('mixedTripletStep'),
    Fitter = cms.string('FlexibleKFFittingSmoother'),
    GeometricInnerState = cms.bool(False),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    Propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    TTRHBuilder = cms.string('WithAngleAndTemplateWithoutProbQ'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("mixedTripletStepTrackCandidates"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(True)
)


process.mkFitEventOfHits = cms.EDProducer("MkFitEventOfHitsProducer",
    beamSpot = cms.InputTag("offlineBeamSpot"),
    mightGet = cms.optional.untracked.vstring,
    pixelHits = cms.InputTag("mkFitSiPixelHits"),
    stripHits = cms.InputTag("mkFitSiStripHits"),
    usePixelQualityDB = cms.bool(True),
    useStripStripQualityDB = cms.bool(True)
)


process.mkFitEventOfHitsPreSplitting = cms.EDProducer("MkFitEventOfHitsProducer",
    beamSpot = cms.InputTag("offlineBeamSpot"),
    mightGet = cms.optional.untracked.vstring,
    pixelHits = cms.InputTag("mkFitSiPixelHitsPreSplitting"),
    stripHits = cms.InputTag("mkFitSiStripHits"),
    usePixelQualityDB = cms.bool(True),
    useStripStripQualityDB = cms.bool(True)
)


process.mkFitSiPixelHits = cms.EDProducer("MkFitSiPixelHitConverter",
    hits = cms.InputTag("siPixelRecHits"),
    mightGet = cms.optional.untracked.vstring,
    ttrhBuilder = cms.ESInputTag("","WithTrackAngle")
)


process.mkFitSiPixelHitsPreSplitting = cms.EDProducer("MkFitSiPixelHitConverter",
    hits = cms.InputTag("siPixelRecHitsPreSplitting"),
    mightGet = cms.optional.untracked.vstring,
    ttrhBuilder = cms.ESInputTag("","WithTrackAngle")
)


process.mkFitSiStripHits = cms.EDProducer("MkFitSiStripHitConverter",
    mightGet = cms.optional.untracked.vstring,
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutLoose')
    ),
    rphiHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
    stereoHits = cms.InputTag("siStripMatchedRecHits","stereoRecHit"),
    ttrhBuilder = cms.ESInputTag("","WithTrackAngle")
)


process.muonSeededSeedsInOut = cms.EDProducer("MuonReSeeder",
    DoPredictionsOnly = cms.bool(False),
    Fitter = cms.string('KFFitterForRefitInsideOut'),
    MTDRecHitBuilder = cms.string('MTDRecHitBuilder'),
    MuonRecHitBuilder = cms.string('MuonRecHitBuilder'),
    Propagator = cms.string('SmartPropagatorAnyRKOpposite'),
    RefitDirection = cms.string('alongMomentum'),
    RefitRPCHits = cms.bool(True),
    Smoother = cms.string('KFSmootherForRefitInsideOut'),
    TrackerRecHitBuilder = cms.string('WithAngleAndTemplate'),
    cut = cms.string('pt > 2'),
    debug = cms.untracked.bool(False),
    insideOut = cms.bool(True),
    layersToKeep = cms.int32(5),
    src = cms.InputTag("earlyMuons")
)


process.muonSeededSeedsInOutAsTracks = cms.EDProducer("FakeTrackProducerFromSeed",
    src = cms.InputTag("muonSeededSeedsInOut")
)


process.muonSeededSeedsOutIn = cms.EDProducer("OutsideInMuonSeeder",
    cut = cms.string('pt > 10 && outerTrack.hitPattern.muonStationsWithValidHits >= 2'),
    debug = cms.untracked.bool(False),
    errorRescaleFactor = cms.double(2.0),
    fromVertex = cms.bool(True),
    hitCollector = cms.string('hitCollectorForOutInMuonSeeds'),
    hitsToTry = cms.int32(3),
    layersToTry = cms.int32(3),
    maxEtaForTOB = cms.double(1.8),
    minEtaForTEC = cms.double(0.7),
    muonPropagator = cms.string('SteppingHelixPropagatorAlong'),
    src = cms.InputTag("earlyMuons"),
    trackerPropagator = cms.string('PropagatorWithMaterial')
)


process.muonSeededSeedsOutInAsTracks = cms.EDProducer("FakeTrackProducerFromSeed",
    src = cms.InputTag("muonSeededSeedsOutIn")
)


process.muonSeededTrackCandidatesInOut = cms.EDProducer("CkfTrackCandidateMaker",
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    RedundantSeedCleaner = cms.string('none'),
    TrajectoryBuilderPSet = cms.PSet(
        refToPSet_ = cms.string('muonSeededTrajectoryBuilderForInOut')
    ),
    TrajectoryCleaner = cms.string('muonSeededTrajectoryCleanerBySharedHits'),
    TransientInitialStateEstimatorParameters = cms.PSet(
        numberMeasurementsForFit = cms.int32(4),
        propagatorAlongTISE = cms.string('PropagatorWithMaterial'),
        propagatorOppositeTISE = cms.string('PropagatorWithMaterialOpposite')
    ),
    cleanTrajectoryAfterInOut = cms.bool(True),
    clustersToSkip = cms.InputTag(""),
    doSeedingRegionRebuilding = cms.bool(True),
    maxNSeeds = cms.uint32(500000),
    maxSeedsBeforeCleaning = cms.uint32(5000),
    mightGet = cms.optional.untracked.vstring,
    numHitsForSeedCleaner = cms.int32(4),
    onlyPixelHitsForSeedCleaner = cms.bool(False),
    phase2clustersToSkip = cms.InputTag(""),
    reverseTrajectories = cms.bool(False),
    src = cms.InputTag("muonSeededSeedsInOut"),
    useHitsSplitting = cms.bool(True)
)


process.muonSeededTrackCandidatesInOutAsTracks = cms.EDProducer("FakeTrackProducerFromCandidate",
    src = cms.InputTag("muonSeededTrackCandidatesInOut")
)


process.muonSeededTrackCandidatesOutIn = cms.EDProducer("CkfTrackCandidateMaker",
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    RedundantSeedCleaner = cms.string('CachingSeedCleanerBySharedInput'),
    TrajectoryBuilderPSet = cms.PSet(
        refToPSet_ = cms.string('muonSeededTrajectoryBuilderForOutIn')
    ),
    TrajectoryCleaner = cms.string('muonSeededTrajectoryCleanerBySharedHits'),
    TransientInitialStateEstimatorParameters = cms.PSet(
        numberMeasurementsForFit = cms.int32(4),
        propagatorAlongTISE = cms.string('PropagatorWithMaterial'),
        propagatorOppositeTISE = cms.string('PropagatorWithMaterialOpposite')
    ),
    cleanTrajectoryAfterInOut = cms.bool(True),
    clustersToSkip = cms.InputTag(""),
    doSeedingRegionRebuilding = cms.bool(True),
    maxNSeeds = cms.uint32(500000),
    maxSeedsBeforeCleaning = cms.uint32(5000),
    mightGet = cms.optional.untracked.vstring,
    numHitsForSeedCleaner = cms.int32(50),
    onlyPixelHitsForSeedCleaner = cms.bool(False),
    phase2clustersToSkip = cms.InputTag(""),
    reverseTrajectories = cms.bool(False),
    src = cms.InputTag("muonSeededSeedsOutIn"),
    useHitsSplitting = cms.bool(True)
)


process.muonSeededTrackCandidatesOutInAsTracks = cms.EDProducer("FakeTrackProducerFromCandidate",
    src = cms.InputTag("muonSeededTrackCandidatesOutIn")
)


process.muonSeededTracksInOut = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('muonSeededStepInOut'),
    Fitter = cms.string('muonSeededFittingSmootherWithOutliersRejectionAndRK'),
    GeometricInnerState = cms.bool(False),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    Propagator = cms.string('RungeKuttaTrackerPropagator'),
    SimpleMagneticField = cms.string(''),
    TTRHBuilder = cms.string('WithAngleAndTemplate'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("muonSeededTrackCandidatesInOut"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(False)
)


process.muonSeededTracksInOutClassifier = cms.EDProducer("TrackCutClassifier",
    beamspot = cms.InputTag("offlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mightGet = cms.optional.untracked.vstring,
    mva = cms.PSet(
        dr_par = cms.PSet(
            d0err = cms.vdouble(0.003, 0.003, 0.003),
            d0err_par = cms.vdouble(0.001, 0.001, 0.001),
            drWPVerr_par = cms.vdouble(3.4028234663852886e+38, 3.4028234663852886e+38, 3.4028234663852886e+38),
            dr_exp = cms.vint32(2147483647, 2147483647, 2147483647),
            dr_par1 = cms.vdouble(3.4028234663852886e+38, 3.4028234663852886e+38, 3.4028234663852886e+38),
            dr_par2 = cms.vdouble(3.4028234663852886e+38, 3.4028234663852886e+38, 3.4028234663852886e+38)
        ),
        dz_par = cms.PSet(
            dzWPVerr_par = cms.vdouble(3.4028234663852886e+38, 3.4028234663852886e+38, 3.4028234663852886e+38),
            dz_exp = cms.vint32(2147483647, 2147483647, 2147483647),
            dz_par1 = cms.vdouble(3.4028234663852886e+38, 3.4028234663852886e+38, 3.4028234663852886e+38),
            dz_par2 = cms.vdouble(3.4028234663852886e+38, 3.4028234663852886e+38, 3.4028234663852886e+38)
        ),
        isHLT = cms.bool(False),
        maxChi2 = cms.vdouble(9999.0, 9999.0, 9999.0),
        maxChi2n = cms.vdouble(10.0, 1.0, 0.4),
        maxDr = cms.vdouble(3.4028234663852886e+38, 3.4028234663852886e+38, 3.4028234663852886e+38),
        maxDz = cms.vdouble(3.4028234663852886e+38, 3.4028234663852886e+38, 3.4028234663852886e+38),
        maxDzWrtBS = cms.vdouble(3.4028234663852886e+38, 24, 15),
        maxLostLayers = cms.vint32(4, 3, 2),
        maxRelPtErr = cms.vdouble(3.4028234663852886e+38, 3.4028234663852886e+38, 3.4028234663852886e+38),
        min3DLayers = cms.vint32(1, 2, 2),
        minHits = cms.vint32(0, 0, 1),
        minHits4pass = cms.vint32(2147483647, 2147483647, 2147483647),
        minLayers = cms.vint32(3, 5, 5),
        minNVtxTrk = cms.int32(2),
        minNdof = cms.vdouble(-1, -1, -1),
        minPixelHits = cms.vint32(0, 0, 0)
    ),
    qualityCuts = cms.vdouble(-0.7, 0.1, 0.7),
    src = cms.InputTag("muonSeededTracksInOut"),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.muonSeededTracksInOutSelector = cms.EDProducer("MultiTrackSelector",
    beamspot = cms.InputTag("offlineBeamSpot"),
    src = cms.InputTag("muonSeededTracksInOut"),
    trackSelectors = cms.VPSet(
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(False),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(10.0),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(0.55, 4.0),
            d0_par2 = cms.vdouble(0.55, 4.0),
            dz_par1 = cms.vdouble(0.65, 4.0),
            dz_par2 = cms.vdouble(0.45, 4.0),
            keepAllTracks = cms.bool(False),
            maxNumberLostLayers = cms.uint32(4),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(7),
            minNumber3DLayers = cms.uint32(0),
            minNumberLayers = cms.uint32(3),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(5),
            nSigmaZ = cms.double(4.0),
            name = cms.string('muonSeededTracksInOutLoose'),
            preFilterName = cms.string(''),
            qualityBit = cms.string('loose'),
            res_par = cms.vdouble(0.003, 0.01),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        ),
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(False),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(1.0),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(0.3, 4.0),
            d0_par2 = cms.vdouble(0.4, 4.0),
            dz_par1 = cms.vdouble(0.35, 4.0),
            dz_par2 = cms.vdouble(0.4, 4.0),
            keepAllTracks = cms.bool(True),
            maxNumberLostLayers = cms.uint32(3),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(10),
            minNumber3DLayers = cms.uint32(2),
            minNumberLayers = cms.uint32(5),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(6),
            nSigmaZ = cms.double(4.0),
            name = cms.string('muonSeededTracksInOutTight'),
            preFilterName = cms.string('muonSeededTracksInOutLoose'),
            qualityBit = cms.string('tight'),
            res_par = cms.vdouble(0.003, 0.01),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        ),
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(False),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(0.4),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(0.3, 4.0),
            d0_par2 = cms.vdouble(0.4, 4.0),
            dz_par1 = cms.vdouble(0.35, 4.0),
            dz_par2 = cms.vdouble(0.4, 4.0),
            keepAllTracks = cms.bool(True),
            maxNumberLostLayers = cms.uint32(2),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(2),
            minNumberLayers = cms.uint32(5),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(7),
            nSigmaZ = cms.double(4.0),
            name = cms.string('muonSeededTracksInOutHighPurity'),
            preFilterName = cms.string('muonSeededTracksInOutTight'),
            qualityBit = cms.string('highPurity'),
            res_par = cms.vdouble(0.003, 0.001),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        )
    ),
    useVertices = cms.bool(True),
    useVtxError = cms.bool(False),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.muonSeededTracksOutIn = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('muonSeededStepOutIn'),
    Fitter = cms.string('muonSeededFittingSmootherWithOutliersRejectionAndRK'),
    GeometricInnerState = cms.bool(False),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    Propagator = cms.string('RungeKuttaTrackerPropagator'),
    SimpleMagneticField = cms.string(''),
    TTRHBuilder = cms.string('WithAngleAndTemplate'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("muonSeededTrackCandidatesOutIn"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(False)
)


process.muonSeededTracksOutInClassifier = cms.EDProducer("TrackCutClassifier",
    beamspot = cms.InputTag("offlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mightGet = cms.optional.untracked.vstring,
    mva = cms.PSet(
        dr_par = cms.PSet(
            d0err = cms.vdouble(0.003, 0.003, 0.003),
            d0err_par = cms.vdouble(0.001, 0.001, 0.001),
            drWPVerr_par = cms.vdouble(3.4028234663852886e+38, 3.4028234663852886e+38, 3.4028234663852886e+38),
            dr_exp = cms.vint32(2147483647, 2147483647, 2147483647),
            dr_par1 = cms.vdouble(3.4028234663852886e+38, 3.4028234663852886e+38, 3.4028234663852886e+38),
            dr_par2 = cms.vdouble(3.4028234663852886e+38, 3.4028234663852886e+38, 3.4028234663852886e+38)
        ),
        dz_par = cms.PSet(
            dzWPVerr_par = cms.vdouble(3.4028234663852886e+38, 3.4028234663852886e+38, 3.4028234663852886e+38),
            dz_exp = cms.vint32(2147483647, 2147483647, 2147483647),
            dz_par1 = cms.vdouble(3.4028234663852886e+38, 3.4028234663852886e+38, 3.4028234663852886e+38),
            dz_par2 = cms.vdouble(3.4028234663852886e+38, 3.4028234663852886e+38, 3.4028234663852886e+38)
        ),
        isHLT = cms.bool(False),
        maxChi2 = cms.vdouble(9999.0, 9999.0, 9999.0),
        maxChi2n = cms.vdouble(10.0, 1.0, 0.4),
        maxDr = cms.vdouble(3.4028234663852886e+38, 3.4028234663852886e+38, 3.4028234663852886e+38),
        maxDz = cms.vdouble(3.4028234663852886e+38, 3.4028234663852886e+38, 3.4028234663852886e+38),
        maxDzWrtBS = cms.vdouble(3.4028234663852886e+38, 24, 15),
        maxLostLayers = cms.vint32(4, 3, 2),
        maxRelPtErr = cms.vdouble(3.4028234663852886e+38, 3.4028234663852886e+38, 3.4028234663852886e+38),
        min3DLayers = cms.vint32(1, 2, 2),
        minHits = cms.vint32(0, 0, 1),
        minHits4pass = cms.vint32(2147483647, 2147483647, 2147483647),
        minLayers = cms.vint32(3, 5, 5),
        minNVtxTrk = cms.int32(2),
        minNdof = cms.vdouble(-1, -1, -1),
        minPixelHits = cms.vint32(0, 0, 0)
    ),
    qualityCuts = cms.vdouble(-0.7, 0.1, 0.7),
    src = cms.InputTag("muonSeededTracksOutIn"),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.muonSeededTracksOutInSelector = cms.EDProducer("MultiTrackSelector",
    beamspot = cms.InputTag("offlineBeamSpot"),
    src = cms.InputTag("muonSeededTracksOutIn"),
    trackSelectors = cms.VPSet(
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(False),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(10.0),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(0.55, 4.0),
            d0_par2 = cms.vdouble(0.55, 4.0),
            dz_par1 = cms.vdouble(0.65, 4.0),
            dz_par2 = cms.vdouble(0.45, 4.0),
            keepAllTracks = cms.bool(False),
            maxNumberLostLayers = cms.uint32(4),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(7),
            minNumber3DLayers = cms.uint32(0),
            minNumberLayers = cms.uint32(3),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(5),
            nSigmaZ = cms.double(4.0),
            name = cms.string('muonSeededTracksOutInLoose'),
            preFilterName = cms.string(''),
            qualityBit = cms.string('loose'),
            res_par = cms.vdouble(0.003, 0.01),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        ),
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(False),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(1.0),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(0.3, 4.0),
            d0_par2 = cms.vdouble(0.4, 4.0),
            dz_par1 = cms.vdouble(0.35, 4.0),
            dz_par2 = cms.vdouble(0.4, 4.0),
            keepAllTracks = cms.bool(True),
            maxNumberLostLayers = cms.uint32(3),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(10),
            minNumber3DLayers = cms.uint32(2),
            minNumberLayers = cms.uint32(5),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(6),
            nSigmaZ = cms.double(4.0),
            name = cms.string('muonSeededTracksOutInTight'),
            preFilterName = cms.string('muonSeededTracksOutInLoose'),
            qualityBit = cms.string('tight'),
            res_par = cms.vdouble(0.003, 0.01),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        ),
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(False),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(0.4),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(0.3, 4.0),
            d0_par2 = cms.vdouble(0.4, 4.0),
            dz_par1 = cms.vdouble(0.35, 4.0),
            dz_par2 = cms.vdouble(0.4, 4.0),
            keepAllTracks = cms.bool(True),
            maxNumberLostLayers = cms.uint32(2),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(2),
            minNumberLayers = cms.uint32(5),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(7),
            nSigmaZ = cms.double(4.0),
            name = cms.string('muonSeededTracksOutInHighPurity'),
            preFilterName = cms.string('muonSeededTracksOutInTight'),
            qualityBit = cms.string('highPurity'),
            res_par = cms.vdouble(0.003, 0.001),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        )
    ),
    useVertices = cms.bool(True),
    useVtxError = cms.bool(False),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.newCombinedSeeds = cms.EDProducer("SeedCombiner",
    seedCollections = cms.VInputTag(
        "initialStepSeeds", "highPtTripletStepSeeds", "mixedTripletStepSeeds", "pixelLessStepSeeds", "tripletElectronSeeds",
        "pixelPairElectronSeeds", "stripPairElectronSeeds", "lowPtTripletStepSeeds", "lowPtQuadStepSeeds", "detachedTripletStepSeeds",
        "detachedQuadStepSeeds", "pixelPairStepSeeds"
    )
)


process.pfRecoTauDiscriminationAgainstElectron = cms.EDProducer("PFRecoTauDiscriminationAgainstElectron",
    ApplyCut_BremCombined = cms.bool(False),
    ApplyCut_BremsRecoveryEOverPLead = cms.bool(False),
    ApplyCut_EOverPLead = cms.bool(False),
    ApplyCut_EcalCrackCut = cms.bool(False),
    ApplyCut_ElectronPreID = cms.bool(False),
    ApplyCut_ElectronPreID_2D = cms.bool(False),
    ApplyCut_EmFraction = cms.bool(False),
    ApplyCut_Hcal3x3OverPLead = cms.bool(False),
    ApplyCut_HcalMaxOverPLead = cms.bool(False),
    ApplyCut_HcalTotOverPLead = cms.bool(False),
    ApplyCut_PFElectronMVA = cms.bool(True),
    BremCombined_Fraction = cms.double(0.99),
    BremCombined_HOP = cms.double(0.1),
    BremCombined_Mass = cms.double(0.55),
    BremCombined_StripSize = cms.double(0.03),
    BremsRecoveryEOverPLead_maxValue = cms.double(1.8),
    BremsRecoveryEOverPLead_minValue = cms.double(0.8),
    EOverPLead_maxValue = cms.double(1.8),
    EOverPLead_minValue = cms.double(0.8),
    ElecPreID0_EOverPLead_maxValue = cms.double(0.95),
    ElecPreID0_HOverPLead_minValue = cms.double(0.05),
    ElecPreID1_EOverPLead_maxValue = cms.double(0.8),
    ElecPreID1_HOverPLead_minValue = cms.double(0.15),
    EmFraction_maxValue = cms.double(0.9),
    Hcal3x3OverPLead_minValue = cms.double(0.1),
    HcalMaxOverPLead_minValue = cms.double(0.1),
    HcalTotOverPLead_minValue = cms.double(0.1),
    PFElectronMVA_maxValue = cms.double(-0.1),
    PFTauProducer = cms.InputTag("pfRecoTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            Producer = cms.InputTag("pfRecoTauDiscriminationByLeadingTrackFinding"),
            cut = cms.double(0.5)
        )
    )
)


process.pfRecoTauDiscriminationAgainstElectronDeadECAL = cms.EDProducer("PFRecoTauDiscriminationAgainstElectronDeadECAL",
    PFTauProducer = cms.InputTag("fixme"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('AND'),
        decayMode = cms.PSet(
            Producer = cms.InputTag("fixme"),
            cut = cms.double(0)
        ),
        leadTrack = cms.PSet(
            Producer = cms.InputTag("fixme"),
            cut = cms.double(0)
        )
    ),
    dR = cms.double(0.08),
    extrapolateToECalEntrance = cms.bool(True),
    mightGet = cms.optional.untracked.vstring,
    minStatus = cms.uint32(12),
    verbosity = cms.int32(0)
)


process.pfRecoTauDiscriminationAgainstElectronMVA6 = cms.EDProducer("PFRecoTauDiscriminationAgainstElectronMVA6",
    PFTauProducer = cms.InputTag("fixme"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('AND'),
        decayMode = cms.PSet(
            Producer = cms.InputTag("fixme"),
            cut = cms.double(0)
        ),
        leadTrack = cms.PSet(
            Producer = cms.InputTag("fixme"),
            cut = cms.double(0)
        )
    ),
    hgcalElectronIDs = cms.VInputTag(),
    isPhase2 = cms.bool(False),
    loadMVAfromDB = cms.bool(True),
    method = cms.string('BDTG'),
    mightGet = cms.optional.untracked.vstring,
    minMVANoEleMatchWOgWOgsfBL = cms.double(0),
    minMVANoEleMatchWOgWOgsfEC = cms.double(0),
    minMVANoEleMatchWgWOgsfBL = cms.double(0),
    minMVANoEleMatchWgWOgsfEC = cms.double(0),
    minMVAWOgWgsfBL = cms.double(0),
    minMVAWOgWgsfEC = cms.double(0),
    minMVAWgWgsfBL = cms.double(0),
    minMVAWgWgsfEC = cms.double(0),
    mvaName_NoEleMatch_wGwoGSF_BL = cms.string('gbr_NoEleMatch_wGwoGSF_BL'),
    mvaName_NoEleMatch_wGwoGSF_EC = cms.string('gbr_NoEleMatch_wGwoGSF_EC'),
    mvaName_NoEleMatch_woGwoGSF_BL = cms.string('gbr_NoEleMatch_woGwoGSF_BL'),
    mvaName_NoEleMatch_woGwoGSF_EC = cms.string('gbr_NoEleMatch_woGwoGSF_EC'),
    mvaName_wGwGSF_BL = cms.string('gbr_wGwGSF_BL'),
    mvaName_wGwGSF_EC = cms.string('gbr_wGwGSF_EC'),
    mvaName_woGwGSF_BL = cms.string('gbr_woGwGSF_BL'),
    mvaName_woGwGSF_EC = cms.string('gbr_woGwGSF_EC'),
    returnMVA = cms.bool(True),
    srcElectrons = cms.InputTag("fixme"),
    usePhiAtEcalEntranceExtrapolation = cms.bool(False),
    verbosity = cms.int32(0),
    vetoEcalCracks = cms.bool(True)
)


process.pfRecoTauDiscriminationAgainstMuon = cms.EDProducer("PFRecoTauDiscriminationAgainstMuon",
    HoPMin = cms.double(0.2),
    PFTauProducer = cms.InputTag("pfRecoTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            Producer = cms.InputTag("pfRecoTauDiscriminationByLeadingTrackFinding"),
            cut = cms.double(0.5)
        )
    ),
    a = cms.double(0.5),
    b = cms.double(0.5),
    c = cms.double(0.0),
    checkNumMatches = cms.bool(False),
    discriminatorOption = cms.string('noSegMatch'),
    maxNumberOfMatches = cms.int32(0)
)


process.pfRecoTauDiscriminationAgainstMuon2Container = cms.EDProducer("PFRecoTauDiscriminationAgainstMuon2Container",
    IDWPdefinitions = cms.VPSet(cms.PSet(
        HoPMin = cms.double(0.2),
        IDname = cms.string('pfRecoTauDiscriminationAgainstMuon2Container'),
        discriminatorOption = cms.string('loose'),
        doCaloMuonVeto = cms.bool(False),
        maxNumberOfHitsLast2Stations = cms.int32(0),
        maxNumberOfMatches = cms.int32(0)
    )),
    IDdefinitions = cms.VPSet(),
    PFTauProducer = cms.InputTag("pfRecoTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            Producer = cms.InputTag("pfRecoTauDiscriminationByLeadingTrackFinding"),
            cut = cms.double(0.5)
        )
    ),
    dRmuonMatch = cms.double(0.3),
    dRmuonMatchLimitedToJetArea = cms.bool(False),
    maskHitsCSC = cms.vint32(0, 0, 0, 0),
    maskHitsDT = cms.vint32(0, 0, 0, 0),
    maskHitsRPC = cms.vint32(0, 0, 0, 0),
    maskMatchesCSC = cms.vint32(1, 0, 0, 0),
    maskMatchesDT = cms.vint32(0, 0, 0, 0),
    maskMatchesRPC = cms.vint32(0, 0, 0, 0),
    mightGet = cms.optional.untracked.vstring,
    minPtMatchedMuon = cms.double(5),
    srcMuons = cms.InputTag("muons"),
    verbosity = cms.int32(0)
)


process.pfRecoTauDiscriminationAgainstMuonMVA = cms.EDProducer("PFRecoTauDiscriminationAgainstMuonMVA",
    PFTauProducer = cms.InputTag("pfTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            Producer = cms.InputTag("pfRecoTauDiscriminationByLeadingTrackFinding"),
            cut = cms.double(0.5)
        )
    ),
    dRmuonMatch = cms.double(0.3),
    inputFileName = cms.FileInPath('RecoTauTag/RecoTau/data/emptyMVAinputFile'),
    loadMVAfromDB = cms.bool(True),
    mvaMin = cms.double(0.0),
    mvaName = cms.string('againstMuonMVA'),
    returnMVA = cms.bool(True),
    srcMuons = cms.InputTag("muons"),
    verbosity = cms.int32(0)
)


process.pfRecoTauDiscriminationByIsolation = cms.EDProducer("PFRecoTauDiscriminationByIsolationContainer",
    IDWPdefinitions = cms.VPSet(),
    IDdefinitions = cms.VPSet(),
    PFTauProducer = cms.InputTag("pfRecoTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        leadTrack = cms.PSet(
            Producer = cms.InputTag("pfRecoTauDiscriminationByLeadingTrackFinding"),
            cut = cms.double(0.5)
        )
    ),
    WeightECALIsolation = cms.double(1),
    applyFootprintCorrection = cms.bool(False),
    applyRhoCorrection = cms.bool(False),
    customOuterCone = cms.double(-1),
    deltaBetaFactor = cms.string('0.38'),
    deltaBetaPUTrackPtCutOverride = cms.bool(False),
    deltaBetaPUTrackPtCutOverride_val = cms.double(-1.5),
    footprintCorrections = cms.VPSet(
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 0')
        ),
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 1 || decayMode() = 2')
        ),
        cms.PSet(
            offset = cms.string('2.7'),
            selection = cms.string('decayMode() = 5')
        ),
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 6')
        ),
        cms.PSet(
            offset = cms.string('max(2.0, 0.22*pt() - 2.0)'),
            selection = cms.string('decayMode() = 10')
        )
    ),
    isoConeSizeForDeltaBeta = cms.double(0.5),
    mightGet = cms.optional.untracked.vstring,
    minTauPtForNoIso = cms.double(-99),
    particleFlowSrc = cms.InputTag("particleFlow"),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.2),
            maxDeltaZToLeadTrack = cms.double(-1.0),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.03),
            minGammaEt = cms.double(1.5),
            minTrackHits = cms.uint32(8),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(1.0),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        leadingTrkOrPFCandOption = cms.string('leadPFCand'),
        primaryVertexSrc = cms.InputTag("offlinePrimaryVertices"),
        pvFindingAlgo = cms.string('closestInDeltaZ'),
        recoverLeadingTrk = cms.bool(False),
        signalQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.4),
            maxDeltaZToLeadTrack = cms.double(-1.0),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(1.0),
            minNeutralHadronEt = cms.double(30.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        vertexTrackFiltering = cms.bool(False),
        vxAssocQualityCuts = cms.PSet(
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(1.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        )
    ),
    rhoConeSize = cms.double(0.5),
    rhoProducer = cms.InputTag("fixedGridRhoFastjetAll"),
    rhoUEOffsetCorrection = cms.double(1),
    verbosity = cms.int32(0),
    vertexSrc = cms.InputTag("offlinePrimaryVertices")
)


process.pfRecoTauDiscriminationByIsolationContainer = cms.EDProducer("PFRecoTauDiscriminationByIsolationContainer",
    IDWPdefinitions = cms.VPSet(),
    IDdefinitions = cms.VPSet(),
    PFTauProducer = cms.InputTag("pfRecoTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        decayMode = cms.PSet(
            Producer = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs"),
            cut = cms.double(0.5)
        ),
        leadTrack = cms.PSet(
            Producer = cms.InputTag("pfRecoTauDiscriminationByLeadingTrackFinding"),
            cut = cms.double(0.5)
        ),
        preIso = cms.PSet(
            Producer = cms.InputTag("hpsPFTauDiscriminationByLooseChargedIsolation"),
            cut = cms.double(0.5)
        )
    ),
    WeightECALIsolation = cms.double(1),
    applyFootprintCorrection = cms.bool(False),
    applyRhoCorrection = cms.bool(False),
    customOuterCone = cms.double(-1),
    deltaBetaFactor = cms.string('0.38'),
    deltaBetaPUTrackPtCutOverride = cms.bool(False),
    deltaBetaPUTrackPtCutOverride_val = cms.double(-1.5),
    footprintCorrections = cms.VPSet(),
    isoConeSizeForDeltaBeta = cms.double(0.5),
    mightGet = cms.optional.untracked.vstring,
    minTauPtForNoIso = cms.double(-99),
    particleFlowSrc = cms.InputTag("particleFlow"),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.2),
            maxDeltaZToLeadTrack = cms.double(-1),
            maxTrackChi2 = cms.double(100),
            maxTransverseImpactParameter = cms.double(0.03),
            minGammaEt = cms.double(1.5),
            minTrackHits = cms.uint32(8),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(1),
            minTrackVertexWeight = cms.double(-1),
            useTracksInsteadOfPFHadrons = cms.optional.bool
        ),
        leadingTrkOrPFCandOption = cms.string('leadPFCand'),
        primaryVertexSrc = cms.InputTag("offlinePrimaryVertices"),
        pvFindingAlgo = cms.string('closestInDeltaZ'),
        recoverLeadingTrk = cms.bool(False),
        signalQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.4),
            maxDeltaZToLeadTrack = cms.double(-1),
            maxTrackChi2 = cms.double(100),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(1),
            minNeutralHadronEt = cms.double(30),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1),
            useTracksInsteadOfPFHadrons = cms.optional.bool
        ),
        vertexTrackFiltering = cms.bool(False),
        vxAssocQualityCuts = cms.PSet(
            maxTrackChi2 = cms.double(100),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(1),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1),
            useTracksInsteadOfPFHadrons = cms.optional.bool
        )
    ),
    rhoConeSize = cms.double(0.5),
    rhoProducer = cms.InputTag("fixedGridRhoFastjetAll"),
    rhoUEOffsetCorrection = cms.double(1),
    verbosity = cms.int32(0),
    vertexSrc = cms.InputTag("offlinePrimaryVertices")
)


process.pfRecoTauDiscriminationByLeadingObjectPtCut = cms.EDProducer("PFRecoTauDiscriminationByLeadingObjectPtCut",
    MinPtLeadingObject = cms.double(5.0),
    PFTauProducer = cms.InputTag("pfRecoTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and')
    ),
    UseOnlyChargedHadrons = cms.bool(False)
)


process.pfRecoTauDiscriminationByLeadingTrackFinding = cms.EDProducer("PFRecoTauDiscriminationByLeadingObjectPtCut",
    MinPtLeadingObject = cms.double(0.0),
    PFTauProducer = cms.InputTag("pfRecoTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and')
    ),
    UseOnlyChargedHadrons = cms.bool(True)
)


process.pfRecoTauDiscriminationByMVAIsolationRun2 = cms.EDProducer("PFRecoTauDiscriminationByMVAIsolationRun2",
    PFTauProducer = cms.InputTag("fixme"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('AND'),
        decayMode = cms.PSet(
            Producer = cms.InputTag("fixme"),
            cut = cms.double(0)
        ),
        leadTrack = cms.PSet(
            Producer = cms.InputTag("fixme"),
            cut = cms.double(0)
        )
    ),
    inputFileName = cms.optional.FileInPath,
    inputIDNameSuffix = cms.string(''),
    loadMVAfromDB = cms.bool(True),
    mightGet = cms.optional.untracked.vstring,
    mvaName = cms.string('tauIdMVAnewDMwLT'),
    mvaOpt = cms.string('newDMwLT'),
    srcBasicTauDiscriminators = cms.InputTag("hpsPFTauBasicDiscriminators"),
    srcTauTransverseImpactParameters = cms.InputTag(""),
    verbosity = cms.int32(0)
)


process.pfRecoTauTagInfoProducer = cms.EDProducer("PFRecoTauTagInfoProducer",
    ChargedHadrCand_AssociationCone = cms.double(0.8),
    ChargedHadrCand_tkPVmaxDZ = cms.double(0.2),
    ChargedHadrCand_tkmaxChi2 = cms.double(100),
    ChargedHadrCand_tkmaxipt = cms.double(0.03),
    ChargedHadrCand_tkminPixelHitsn = cms.int32(0),
    ChargedHadrCand_tkminPt = cms.double(0.5),
    ChargedHadrCand_tkminTrackerHitsn = cms.int32(3),
    GammaCand_EcalclusMinEt = cms.double(1),
    NeutrHadrCand_HcalclusMinEt = cms.double(1),
    PFCandidateProducer = cms.InputTag("particleFlow"),
    PFJetTracksAssociatorProducer = cms.InputTag("ak4PFJetTracksAssociatorAtVertex"),
    PVProducer = cms.InputTag("offlinePrimaryVertices"),
    UsePVconstraint = cms.bool(True),
    mightGet = cms.optional.untracked.vstring,
    smearedPVsigmaX = cms.double(0.0015),
    smearedPVsigmaY = cms.double(0.0015),
    smearedPVsigmaZ = cms.double(0.005),
    tkPVmaxDZ = cms.double(0.2),
    tkmaxChi2 = cms.double(100),
    tkmaxipt = cms.double(0.03),
    tkminPixelHitsn = cms.int32(0),
    tkminPt = cms.double(0.5),
    tkminTrackerHitsn = cms.int32(3)
)


process.photonConvTrajSeedFromQuadruplets = cms.EDProducer("PhotonConversionTrajectorySeedProducerFromQuadruplets",
    ClusterCheckPSet = cms.PSet(
        ClusterCollectionLabel = cms.InputTag("siStripClusters"),
        MaxNumberOfCosmicClusters = cms.uint32(50000),
        MaxNumberOfPixelClusters = cms.uint32(10000),
        PixelClusterCollectionLabel = cms.InputTag("siPixelClusters"),
        doClusterCheck = cms.bool(True)
    ),
    DoxcheckSeedCandidates = cms.bool(False),
    OrderedHitsFactoryPSet = cms.PSet(
        SeedingLayers = cms.InputTag("conv2LayerPairs"),
        maxElement = cms.uint32(900000)
    ),
    QuadCutPSet = cms.PSet(
        Cut_BeamPipeRadius = cms.double(3.0),
        Cut_DeltaRho = cms.double(12.0),
        Cut_maxLegPt = cms.double(10.0),
        Cut_minLegPt = cms.double(0.6),
        Cut_zCA = cms.double(100),
        apply_Arbitration = cms.bool(True),
        apply_ClusterShapeFilter = cms.bool(True),
        apply_DeltaPhiCuts = cms.bool(True),
        apply_zCACut = cms.bool(False),
        rejectAllQuads = cms.bool(False)
    ),
    RegionFactoryPSet = cms.PSet(
        ComponentName = cms.string('GlobalRegionProducerFromBeamSpot'),
        RegionPSet = cms.PSet(
            beamSpot = cms.InputTag("offlineBeamSpot"),
            originHalfLength = cms.double(12.0),
            originRadius = cms.double(3.0),
            precise = cms.bool(True),
            ptMin = cms.double(0.2)
        )
    ),
    SeedComparitorPSet = cms.PSet(
        ClusterShapeCacheSrc = cms.InputTag("siPixelClusterShapeCache"),
        ClusterShapeHitFilterName = cms.string('ClusterShapeHitFilter'),
        ComponentName = cms.string('PixelClusterShapeSeedComparitor'),
        FilterAtHelixStage = cms.bool(True),
        FilterPixelHits = cms.bool(False),
        FilterStripHits = cms.bool(True)
    ),
    SeedCreatorPSet = cms.PSet(
        ComponentName = cms.string('SeedForPhotonConversionFromQuadruplets'),
        SeedMomentumForBOFF = cms.double(5.0),
        propagator = cms.string('PropagatorWithMaterial')
    ),
    TrackRefitter = cms.InputTag("generalTracks"),
    beamSpotInputTag = cms.InputTag("offlineBeamSpot"),
    newSeedCandidates = cms.string('conv2SeedCandidates'),
    primaryVerticesTag = cms.InputTag("pixelVertices"),
    xcheckSeedCandidates = cms.string('xcheckSeedCandidates')
)


process.photonConvTrajSeedFromSingleLeg = cms.EDProducer("PhotonConversionTrajectorySeedProducerFromSingleLeg",
    ClusterCheckPSet = cms.PSet(
        ClusterCollectionLabel = cms.InputTag("siStripClusters"),
        MaxNumberOfCosmicClusters = cms.uint32(400000),
        MaxNumberOfPixelClusters = cms.uint32(40000),
        PixelClusterCollectionLabel = cms.InputTag("siPixelClusters"),
        cut = cms.string('strip < 400000 && pixel < 40000 && (strip < 50000 + 10*pixel) && (pixel < 5000 + 0.1*strip)'),
        doClusterCheck = cms.bool(False)
    ),
    DoxcheckSeedCandidates = cms.bool(False),
    OrderedHitsFactoryPSet = cms.PSet(
        SeedingLayers = cms.InputTag("convLayerPairs"),
        maxElement = cms.uint32(40000),
        maxHitPairsPerTrackAndGenerator = cms.uint32(10)
    ),
    RegionFactoryPSet = cms.PSet(
        ComponentName = cms.string('GlobalRegionProducerFromBeamSpot'),
        RegionPSet = cms.PSet(
            beamSpot = cms.InputTag("offlineBeamSpot"),
            originHalfLength = cms.double(12.0),
            originRadius = cms.double(3.0),
            precise = cms.bool(True),
            ptMin = cms.double(0.2)
        )
    ),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('none')
    ),
    SeedCreatorPSet = cms.PSet(
        ComponentName = cms.string('SeedForPhotonConversion1Leg'),
        SeedMomentumForBOFF = cms.double(5.0),
        TTRHBuilder = cms.string('WithTrackAngle'),
        propagator = cms.string('PropagatorWithMaterial')
    ),
    TrackRefitter = cms.InputTag("generalTracks"),
    applyTkVtxConstraint = cms.bool(True),
    beamSpotInputTag = cms.InputTag("offlineBeamSpot"),
    maxDZSigmas = cms.double(10.0),
    maxNumSelVtx = cms.uint32(2),
    newSeedCandidates = cms.string('convSeedCandidates'),
    primaryVerticesTag = cms.InputTag("firstStepPrimaryVertices"),
    vtxMinDoF = cms.double(4),
    xcheckSeedCandidates = cms.string('xcheckSeedCandidates')
)


process.pixelLessLayerPairs4PixelLessTracking = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(

    ),
    FPix = cms.PSet(

    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        maxRing = cms.int32(2),
        minRing = cms.int32(1),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHitUnmatched"),
        stereoRecHits = cms.InputTag("siStripMatchedRecHits","stereoRecHitUnmatched"),
        useRingSlector = cms.bool(True),
        useSimpleRphiHitsCleaner = cms.bool(False)
    ),
    TIB = cms.PSet(

    ),
    TIB1 = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHitUnmatched"),
        stereoRecHits = cms.InputTag("siStripMatchedRecHits","stereoRecHitUnmatched"),
        useSimpleRphiHitsCleaner = cms.bool(False)
    ),
    TIB2 = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHitUnmatched"),
        stereoRecHits = cms.InputTag("siStripMatchedRecHits","stereoRecHitUnmatched"),
        useSimpleRphiHitsCleaner = cms.bool(False)
    ),
    TIB3 = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHitUnmatched"),
        stereoRecHits = cms.InputTag("siStripMatchedRecHits","stereoRecHitUnmatched"),
        useSimpleRphiHitsCleaner = cms.bool(False)
    ),
    TID = cms.PSet(

    ),
    TID1 = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        maxRing = cms.int32(3),
        minRing = cms.int32(1),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHitUnmatched"),
        stereoRecHits = cms.InputTag("siStripMatchedRecHits","stereoRecHitUnmatched"),
        useRingSlector = cms.bool(True),
        useSimpleRphiHitsCleaner = cms.bool(False)
    ),
    TID2 = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        maxRing = cms.int32(3),
        minRing = cms.int32(1),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHitUnmatched"),
        stereoRecHits = cms.InputTag("siStripMatchedRecHits","stereoRecHitUnmatched"),
        useRingSlector = cms.bool(True),
        useSimpleRphiHitsCleaner = cms.bool(False)
    ),
    TID3 = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        maxRing = cms.int32(2),
        minRing = cms.int32(1),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHitUnmatched"),
        stereoRecHits = cms.InputTag("siStripMatchedRecHits","stereoRecHitUnmatched"),
        useRingSlector = cms.bool(True),
        useSimpleRphiHitsCleaner = cms.bool(False)
    ),
    TOB = cms.PSet(

    ),
    layerList = cms.vstring(
        'TIB1+TIB2',
        'TIB1+TIB3',
        'TIB2+TIB3',
        'TIB1+TID1_pos',
        'TIB1+TID1_neg',
        'TIB2+TID1_pos',
        'TIB2+TID1_neg',
        'TIB1+TID2_pos',
        'TIB1+TID2_neg',
        'TID1_pos+TID2_pos',
        'TID2_pos+TID3_pos',
        'TID3_pos+TEC2_pos',
        'TEC1_pos+TEC2_pos',
        'TEC2_pos+TEC3_pos',
        'TID1_neg+TID2_neg',
        'TID2_neg+TID3_neg',
        'TID3_neg+TEC2_neg',
        'TEC1_neg+TEC2_neg',
        'TEC2_neg+TEC3_neg'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.pixelLessStep = cms.EDProducer("TrackTfClassifier",
    beamspot = cms.InputTag("offlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mightGet = cms.optional.untracked.vstring,
    mva = cms.PSet(
        batchSize = cms.int32(16),
        tfDnnLabel = cms.string('trackSelectionTf_CKF')
    ),
    qualityCuts = cms.vdouble(-0.82, -0.61, -0.16),
    src = cms.InputTag("pixelLessStepTracks"),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.pixelLessStepClassifier1 = cms.EDProducer("TrackMVAClassifierPrompt",
    beamspot = cms.InputTag("offlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mightGet = cms.optional.untracked.vstring,
    mva = cms.PSet(
        GBRForestFileName = cms.string(''),
        GBRForestLabel = cms.string('MVASelectorIter5_13TeV')
    ),
    qualityCuts = cms.vdouble(-0.4, 0.0, 0.4),
    src = cms.InputTag("pixelLessStepTracks"),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.pixelLessStepClassifier2 = cms.EDProducer("TrackMVAClassifierPrompt",
    beamspot = cms.InputTag("offlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mightGet = cms.optional.untracked.vstring,
    mva = cms.PSet(
        GBRForestFileName = cms.string(''),
        GBRForestLabel = cms.string('MVASelectorIter0_13TeV')
    ),
    qualityCuts = cms.vdouble(-0.0, 0.0, 0.0),
    src = cms.InputTag("pixelLessStepTracks"),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.pixelLessStepClusters = cms.EDProducer("TrackClusterRemover",
    TrackQuality = cms.string('highPurity'),
    maxChi2 = cms.double(9.0),
    mightGet = cms.optional.untracked.vstring,
    minNumberOfLayersWithMeasBeforeFiltering = cms.int32(0),
    oldClusterRemovalInfo = cms.InputTag("mixedTripletStepClusters"),
    overrideTrkQuals = cms.InputTag(""),
    pixelClusters = cms.InputTag("siPixelClusters"),
    stripClusters = cms.InputTag("siStripClusters"),
    trackClassifier = cms.InputTag("mixedTripletStep","QualityMasks"),
    trajectories = cms.InputTag("mixedTripletStepTracks")
)


process.pixelLessStepHitDoublets = cms.EDProducer("HitPairEDProducer",
    clusterCheck = cms.InputTag("trackerClusterCheck"),
    layerPairs = cms.vuint32(0),
    maxElement = cms.uint32(50000000),
    maxElementTotal = cms.uint32(50000000),
    mightGet = cms.optional.untracked.vstring,
    produceIntermediateHitDoublets = cms.bool(True),
    produceSeedingHitSets = cms.bool(False),
    seedingLayers = cms.InputTag("pixelLessStepSeedLayers"),
    trackingRegions = cms.InputTag("pixelLessStepTrackingRegions"),
    trackingRegionsSeedingLayers = cms.InputTag("")
)


process.pixelLessStepHitTriplets = cms.EDProducer("MultiHitFromChi2EDProducer",
    ClusterShapeHitFilterName = cms.string('ClusterShapeHitFilter'),
    TTRHBuilder = cms.string('WithTrackAngle'),
    chi2VsPtCut = cms.bool(True),
    chi2_cuts = cms.vdouble(3, 4, 5, 5),
    detIdsToDebug = cms.vint32(0, 0, 0),
    doublets = cms.InputTag("pixelLessStepHitDoublets"),
    extraHitRPhitolerance = cms.double(0),
    extraHitRZtolerance = cms.double(0),
    extraPhiKDBox = cms.double(0.005),
    extraRKDBox = cms.double(0.2),
    extraZKDBox = cms.double(0.2),
    fnSigmaRZ = cms.double(2),
    maxChi2 = cms.double(5),
    maxElement = cms.uint32(1000000),
    mightGet = cms.optional.untracked.vstring,
    phiPreFiltering = cms.double(0.3),
    pt_interv = cms.vdouble(0.4, 0.7, 1, 2),
    refitHits = cms.bool(True),
    useFixedPreFiltering = cms.bool(False)
)


process.pixelLessStepMasks = cms.EDProducer("FastTrackerRecHitMaskProducer",
    TrackQuality = cms.string('highPurity'),
    minNumberOfLayersWithMeasBeforeFiltering = cms.int32(0),
    minNumberOfLayersWithMeasBeforeFiltering_ = cms.int32(0),
    oldHitRemovalInfo = cms.InputTag("mixedTripletStepMasks"),
    recHits = cms.InputTag("fastTrackerRecHits"),
    trackClassifier = cms.InputTag("mixedTripletStep","QualityMasks"),
    trackQuality = cms.string('hightPurity'),
    trajectories = cms.InputTag("mixedTripletStepTracks")
)


process.pixelLessStepSeedClusterMask = cms.EDProducer("SeedClusterRemover",
    Common = cms.PSet(
        maxChi2 = cms.double(9.0)
    ),
    oldClusterRemovalInfo = cms.InputTag("mixedTripletStepSeedClusterMask"),
    pixelClusters = cms.InputTag("siPixelClusters"),
    stripClusters = cms.InputTag("siStripClusters"),
    trajectories = cms.InputTag("pixelLessStepSeeds")
)


process.pixelLessStepSeedLayers = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(

    ),
    FPix = cms.PSet(

    ),
    MTEC = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTight')
        ),
        maxRing = cms.int32(3),
        minRing = cms.int32(3),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
        skipClusters = cms.InputTag("pixelLessStepClusters"),
        useRingSlector = cms.bool(True)
    ),
    MTIB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTight')
        ),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
        skipClusters = cms.InputTag("pixelLessStepClusters")
    ),
    MTID = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTight')
        ),
        maxRing = cms.int32(3),
        minRing = cms.int32(3),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
        skipClusters = cms.InputTag("pixelLessStepClusters"),
        useRingSlector = cms.bool(True)
    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTight')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        maxRing = cms.int32(2),
        minRing = cms.int32(1),
        skipClusters = cms.InputTag("pixelLessStepClusters"),
        useRingSlector = cms.bool(True)
    ),
    TIB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTight')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        skipClusters = cms.InputTag("pixelLessStepClusters")
    ),
    TID = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTight')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        maxRing = cms.int32(2),
        minRing = cms.int32(1),
        skipClusters = cms.InputTag("pixelLessStepClusters"),
        useRingSlector = cms.bool(True)
    ),
    TOB = cms.PSet(

    ),
    layerList = cms.vstring(
        'TIB1+TIB2+MTIB3',
        'TIB1+TIB2+MTIB4',
        'TIB1+TIB2+MTID1_pos',
        'TIB1+TIB2+MTID1_neg',
        'TID1_pos+TID2_pos+TID3_pos',
        'TID1_neg+TID2_neg+TID3_neg',
        'TID1_pos+TID2_pos+MTID3_pos',
        'TID1_neg+TID2_neg+MTID3_neg',
        'TID1_pos+TID2_pos+MTEC1_pos',
        'TID1_neg+TID2_neg+MTEC1_neg',
        'TID2_pos+TID3_pos+TEC1_pos',
        'TID2_neg+TID3_neg+TEC1_neg',
        'TID2_pos+TID3_pos+MTEC1_pos',
        'TID2_neg+TID3_neg+MTEC1_neg',
        'TEC1_pos+TEC2_pos+TEC3_pos',
        'TEC1_neg+TEC2_neg+TEC3_neg',
        'TEC1_pos+TEC2_pos+MTEC3_pos',
        'TEC1_neg+TEC2_neg+MTEC3_neg',
        'TEC1_pos+TEC2_pos+TEC4_pos',
        'TEC1_neg+TEC2_neg+TEC4_neg',
        'TEC1_pos+TEC2_pos+MTEC4_pos',
        'TEC1_neg+TEC2_neg+MTEC4_neg',
        'TEC2_pos+TEC3_pos+TEC4_pos',
        'TEC2_neg+TEC3_neg+TEC4_neg',
        'TEC2_pos+TEC3_pos+MTEC4_pos',
        'TEC2_neg+TEC3_neg+MTEC4_neg',
        'TEC2_pos+TEC3_pos+TEC5_pos',
        'TEC2_neg+TEC3_neg+TEC5_neg',
        'TEC2_pos+TEC3_pos+TEC6_pos',
        'TEC2_neg+TEC3_neg+TEC6_neg',
        'TEC3_pos+TEC4_pos+TEC5_pos',
        'TEC3_neg+TEC4_neg+TEC5_neg',
        'TEC3_pos+TEC4_pos+MTEC5_pos',
        'TEC3_neg+TEC4_neg+MTEC5_neg',
        'TEC3_pos+TEC5_pos+TEC6_pos',
        'TEC3_neg+TEC5_neg+TEC6_neg',
        'TEC4_pos+TEC5_pos+TEC6_pos',
        'TEC4_neg+TEC5_neg+TEC6_neg'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.pixelLessStepSeeds = cms.EDProducer("SeedCreatorFromRegionConsecutiveHitsTripletOnlyEDProducer",
    MinOneOverPtError = cms.double(1),
    OriginTransverseErrorMultiplier = cms.double(1),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('CombinedSeedComparitor'),
        comparitors = cms.VPSet(
            cms.PSet(
                ClusterShapeCacheSrc = cms.InputTag("siPixelClusterShapeCache"),
                ClusterShapeHitFilterName = cms.string('pixelLessStepClusterShapeHitFilter'),
                ComponentName = cms.string('PixelClusterShapeSeedComparitor'),
                FilterAtHelixStage = cms.bool(True),
                FilterPixelHits = cms.bool(False),
                FilterStripHits = cms.bool(True)
            ),
            cms.PSet(
                ComponentName = cms.string('StripSubClusterShapeSeedFilter'),
                FilterAtHelixStage = cms.bool(False),
                label = cms.untracked.string('Seeds'),
                layerMask = cms.PSet(

                ),
                maxNSat = cms.uint32(3),
                maxTrimmedSizeDiffNeg = cms.double(1.0),
                maxTrimmedSizeDiffPos = cms.double(0.7),
                seedCutMIPs = cms.double(0.35),
                seedCutSN = cms.double(7.0),
                subclusterCutMIPs = cms.double(0.45),
                subclusterCutSN = cms.double(12.0),
                subclusterWindow = cms.double(0.7),
                trimMaxADC = cms.double(30.0),
                trimMaxFracNeigh = cms.double(0.25),
                trimMaxFracTotal = cms.double(0.15)
            )
        ),
        mode = cms.string('and')
    ),
    SeedMomentumForBOFF = cms.double(5),
    TTRHBuilder = cms.string('WithTrackAngle'),
    forceKinematicWithRegionDirection = cms.bool(False),
    magneticField = cms.string('ParabolicMf'),
    mightGet = cms.optional.untracked.vstring,
    propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    seedingHitSets = cms.InputTag("pixelLessStepHitTriplets")
)


process.pixelLessStepSeeds_vectorHits = cms.EDProducer("SeedingOTEDProducer",
    beamSpotLabel = cms.InputTag("offlineBeamSpot"),
    mightGet = cms.optional.untracked.vstring,
    src = cms.InputTag("siPhase2VectorHits","accepted"),
    trackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    updator = cms.string('KFUpdator')
)


process.pixelLessStepSelector = cms.EDProducer("MultiTrackSelector",
    GBRForestLabel = cms.string('MVASelectorIter5'),
    beamspot = cms.InputTag("offlineBeamSpot"),
    src = cms.InputTag("pixelLessStepTracks"),
    trackSelectors = cms.VPSet(
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(True),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(0.5),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(1.3, 4.0),
            d0_par2 = cms.vdouble(1.3, 4.0),
            dz_par1 = cms.vdouble(1.3, 4.0),
            dz_par2 = cms.vdouble(1.3, 4.0),
            keepAllTracks = cms.bool(False),
            maxNumberLostLayers = cms.uint32(1),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(3),
            minNumberLayers = cms.uint32(4),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(0),
            nSigmaZ = cms.double(4.0),
            name = cms.string('pixelLessStepLoose'),
            preFilterName = cms.string(''),
            qualityBit = cms.string('loose'),
            res_par = cms.vdouble(0.003, 0.001),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        ),
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(True),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(0.35),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(1.1, 4.0),
            d0_par2 = cms.vdouble(1.1, 4.0),
            dz_par1 = cms.vdouble(1.1, 4.0),
            dz_par2 = cms.vdouble(1.1, 4.0),
            keepAllTracks = cms.bool(True),
            maxNumberLostLayers = cms.uint32(0),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(3),
            minNumberLayers = cms.uint32(4),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(0),
            nSigmaZ = cms.double(4.0),
            name = cms.string('pixelLessStepTight'),
            preFilterName = cms.string('pixelLessStepLoose'),
            qualityBit = cms.string('tight'),
            res_par = cms.vdouble(0.003, 0.001),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        ),
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(True),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(0.2),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(0.9, 4.0),
            d0_par2 = cms.vdouble(0.9, 4.0),
            dz_par1 = cms.vdouble(0.9, 4.0),
            dz_par2 = cms.vdouble(0.9, 4.0),
            keepAllTracks = cms.bool(True),
            maxNumberLostLayers = cms.uint32(0),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(3),
            minNumberLayers = cms.uint32(4),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(0),
            nSigmaZ = cms.double(4.0),
            name = cms.string('QualityMasks'),
            preFilterName = cms.string('pixelLessStepTight'),
            qualityBit = cms.string('highPurity'),
            res_par = cms.vdouble(0.003, 0.001),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        )
    ),
    useAnyMVA = cms.bool(False),
    useVertices = cms.bool(True),
    useVtxError = cms.bool(False),
    vertices = cms.InputTag("pixelVertices")
)


process.pixelLessStepTrackCandidates = cms.EDProducer("CkfTrackCandidateMaker",
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    RedundantSeedCleaner = cms.string('CachingSeedCleanerBySharedInput'),
    TrajectoryBuilderPSet = cms.PSet(
        refToPSet_ = cms.string('pixelLessStepTrajectoryBuilder')
    ),
    TrajectoryCleaner = cms.string('pixelLessStepTrajectoryCleanerBySharedHits'),
    TransientInitialStateEstimatorParameters = cms.PSet(
        numberMeasurementsForFit = cms.int32(4),
        propagatorAlongTISE = cms.string('PropagatorWithMaterial'),
        propagatorOppositeTISE = cms.string('PropagatorWithMaterialOpposite')
    ),
    cleanTrajectoryAfterInOut = cms.bool(True),
    clustersToSkip = cms.InputTag("pixelLessStepClusters"),
    doSeedingRegionRebuilding = cms.bool(True),
    maxNSeeds = cms.uint32(500000),
    maxSeedsBeforeCleaning = cms.uint32(5000),
    mightGet = cms.optional.untracked.vstring,
    numHitsForSeedCleaner = cms.int32(50),
    onlyPixelHitsForSeedCleaner = cms.bool(False),
    phase2clustersToSkip = cms.InputTag(""),
    reverseTrajectories = cms.bool(False),
    src = cms.InputTag("pixelLessStepSeeds"),
    useHitsSplitting = cms.bool(True)
)


process.pixelLessStepTrackCandidatesMkFit = cms.EDProducer("MkFitProducer",
    backwardFitInCMSSW = cms.bool(False),
    buildingRoutine = cms.string('cloneEngine'),
    clustersToSkip = cms.InputTag("pixelLessStepClusters"),
    config = cms.ESInputTag("","pixelLessStepTrackCandidatesMkFitConfig"),
    eventOfHits = cms.InputTag("mkFitEventOfHits"),
    limitConcurrency = cms.untracked.bool(False),
    mightGet = cms.optional.untracked.vstring,
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutLoose')
    ),
    mkFitSilent = cms.untracked.bool(True),
    pixelHits = cms.InputTag("mkFitSiPixelHits"),
    removeDuplicates = cms.bool(True),
    seedCleaning = cms.bool(True),
    seeds = cms.InputTag("pixelLessStepTrackCandidatesMkFitSeeds"),
    stripHits = cms.InputTag("mkFitSiStripHits")
)


process.pixelLessStepTrackCandidatesMkFitSeeds = cms.EDProducer("MkFitSeedConverter",
    maxNSeeds = cms.uint32(500000),
    mightGet = cms.optional.untracked.vstring,
    seeds = cms.InputTag("pixelLessStepSeeds"),
    ttrhBuilder = cms.ESInputTag("","WithTrackAngle")
)


process.pixelLessStepTrackingRegions = cms.EDProducer("GlobalTrackingRegionFromBeamSpotEDProducer",
    RegionPSet = cms.PSet(
        beamSpot = cms.InputTag("offlineBeamSpot"),
        nSigmaZ = cms.double(0),
        originHalfLength = cms.double(12.0),
        originRadius = cms.double(1.0),
        precise = cms.bool(True),
        ptMin = cms.double(0.4),
        useMultipleScattering = cms.bool(False)
    ),
    mightGet = cms.optional.untracked.vstring
)


process.pixelLessStepTracks = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('pixelLessStep'),
    Fitter = cms.string('FlexibleKFFittingSmoother'),
    GeometricInnerState = cms.bool(False),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    Propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    TTRHBuilder = cms.string('WithAngleAndTemplateWithoutProbQ'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("pixelLessStepTrackCandidates"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(True)
)


process.pixelPairElectronHitDoublets = cms.EDProducer("HitPairEDProducer",
    clusterCheck = cms.InputTag("trackerClusterCheck"),
    layerPairs = cms.vuint32(0),
    maxElement = cms.uint32(1000000),
    maxElementTotal = cms.uint32(12000000),
    mightGet = cms.optional.untracked.vstring,
    produceIntermediateHitDoublets = cms.bool(False),
    produceSeedingHitSets = cms.bool(True),
    seedingLayers = cms.InputTag("pixelPairElectronSeedLayers"),
    trackingRegions = cms.InputTag("pixelPairElectronTrackingRegions"),
    trackingRegionsSeedingLayers = cms.InputTag("")
)


process.pixelPairElectronSeedLayers = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(
        HitProducer = cms.string('siPixelRecHits'),
        TTRHBuilder = cms.string('WithTrackAngle'),
        skipClusters = cms.InputTag("tripletElectronClusterMask")
    ),
    FPix = cms.PSet(
        HitProducer = cms.string('siPixelRecHits'),
        TTRHBuilder = cms.string('WithTrackAngle'),
        skipClusters = cms.InputTag("tripletElectronClusterMask")
    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(

    ),
    TIB = cms.PSet(

    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(

    ),
    layerList = cms.vstring(
        'BPix1+BPix2',
        'BPix1+BPix3',
        'BPix1+BPix4',
        'BPix2+BPix3',
        'BPix2+BPix4',
        'BPix3+BPix4',
        'BPix1+FPix1_pos',
        'BPix1+FPix1_neg',
        'BPix1+FPix2_pos',
        'BPix1+FPix2_neg',
        'BPix2+FPix1_pos',
        'BPix2+FPix1_neg',
        'FPix1_pos+FPix2_pos',
        'FPix1_neg+FPix2_neg',
        'FPix1_pos+FPix3_pos',
        'FPix1_neg+FPix3_neg',
        'FPix2_pos+FPix3_pos',
        'FPix2_neg+FPix3_neg'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.pixelPairElectronSeeds = cms.EDProducer("SeedCreatorFromRegionConsecutiveHitsEDProducer",
    MinOneOverPtError = cms.double(1),
    OriginTransverseErrorMultiplier = cms.double(1),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('none')
    ),
    SeedMomentumForBOFF = cms.double(5),
    TTRHBuilder = cms.string('WithTrackAngle'),
    forceKinematicWithRegionDirection = cms.bool(False),
    magneticField = cms.string('ParabolicMf'),
    mightGet = cms.optional.untracked.vstring,
    propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    seedingHitSets = cms.InputTag("pixelPairElectronHitDoublets")
)


process.pixelPairElectronTrackingRegions = cms.EDProducer("GlobalTrackingRegionWithVerticesEDProducer",
    RegionPSet = cms.PSet(
        VertexCollection = cms.InputTag("firstStepPrimaryVertices"),
        beamSpot = cms.InputTag("offlineBeamSpot"),
        fixedError = cms.double(0.03),
        halfLengthScaling4BigEvts = cms.bool(False),
        maxNVertices = cms.int32(-1),
        maxPtMin = cms.double(1000),
        minHalfLength = cms.double(0),
        minOriginR = cms.double(0),
        nSigmaZ = cms.double(4),
        originRScaling4BigEvts = cms.bool(False),
        originRadius = cms.double(0.015),
        pixelClustersForScaling = cms.InputTag("siPixelClusters"),
        precise = cms.bool(True),
        ptMin = cms.double(1.0),
        ptMinScaling4BigEvts = cms.bool(False),
        scalingEndNPix = cms.double(1),
        scalingStartNPix = cms.double(0),
        sigmaZVertex = cms.double(3),
        useFakeVertices = cms.bool(False),
        useFixedError = cms.bool(True),
        useFoundVertices = cms.bool(True),
        useMultipleScattering = cms.bool(False)
    ),
    mightGet = cms.optional.untracked.vstring
)


process.pixelPairStep = cms.EDProducer("TrackTfClassifier",
    beamspot = cms.InputTag("offlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mightGet = cms.optional.untracked.vstring,
    mva = cms.PSet(
        batchSize = cms.int32(16),
        tfDnnLabel = cms.string('trackSelectionTf')
    ),
    qualityCuts = cms.vdouble(-0.31, -0.13, 0.13),
    src = cms.InputTag("pixelPairStepTracks"),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.pixelPairStepClusters = cms.EDProducer("TrackClusterRemover",
    TrackQuality = cms.string('highPurity'),
    maxChi2 = cms.double(9.0),
    mightGet = cms.optional.untracked.vstring,
    minNumberOfLayersWithMeasBeforeFiltering = cms.int32(0),
    oldClusterRemovalInfo = cms.InputTag("detachedTripletStepClusters"),
    overrideTrkQuals = cms.InputTag(""),
    pixelClusters = cms.InputTag("siPixelClusters"),
    stripClusters = cms.InputTag("siStripClusters"),
    trackClassifier = cms.InputTag("detachedTripletStep","QualityMasks"),
    trajectories = cms.InputTag("detachedTripletStepTracks")
)


process.pixelPairStepHitDoublets = cms.EDProducer("HitPairEDProducer",
    clusterCheck = cms.InputTag("trackerClusterCheck"),
    layerPairs = cms.vuint32(0),
    maxElement = cms.uint32(1000000),
    maxElementTotal = cms.uint32(12000000),
    mightGet = cms.optional.untracked.vstring,
    produceIntermediateHitDoublets = cms.bool(False),
    produceSeedingHitSets = cms.bool(True),
    seedingLayers = cms.InputTag("pixelPairStepSeedLayers"),
    trackingRegions = cms.InputTag("pixelPairStepTrackingRegions"),
    trackingRegionsSeedingLayers = cms.InputTag("")
)


process.pixelPairStepHitDoubletsB = cms.EDProducer("HitPairEDProducer",
    clusterCheck = cms.InputTag("trackerClusterCheck"),
    layerPairs = cms.vuint32(0),
    maxElement = cms.uint32(1000000),
    maxElementTotal = cms.uint32(12000000),
    mightGet = cms.optional.untracked.vstring,
    produceIntermediateHitDoublets = cms.bool(False),
    produceSeedingHitSets = cms.bool(True),
    seedingLayers = cms.InputTag(""),
    trackingRegions = cms.InputTag(""),
    trackingRegionsSeedingLayers = cms.InputTag("pixelPairStepTrackingRegionsSeedLayersB")
)


process.pixelPairStepMasks = cms.EDProducer("FastTrackerRecHitMaskProducer",
    TrackQuality = cms.string('highPurity'),
    minNumberOfLayersWithMeasBeforeFiltering = cms.int32(0),
    minNumberOfLayersWithMeasBeforeFiltering_ = cms.int32(0),
    oldHitRemovalInfo = cms.InputTag("detachedTripletStepMasks"),
    recHits = cms.InputTag("fastTrackerRecHits"),
    trackClassifier = cms.InputTag("detachedTripletStep","QualityMasks"),
    trackQuality = cms.string('hightPurity'),
    trajectories = cms.InputTag("detachedTripletStepTracks")
)


process.pixelPairStepSeedClusterMask = cms.EDProducer("SeedClusterRemover",
    Common = cms.PSet(
        maxChi2 = cms.double(9.0)
    ),
    oldClusterRemovalInfo = cms.InputTag("initialStepSeedClusterMask"),
    pixelClusters = cms.InputTag("siPixelClusters"),
    stripClusters = cms.InputTag("siStripClusters"),
    trajectories = cms.InputTag("pixelPairStepSeeds")
)


process.pixelPairStepSeedLayers = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(
        HitProducer = cms.string('siPixelRecHits'),
        TTRHBuilder = cms.string('WithTrackAngle'),
        skipClusters = cms.InputTag("pixelPairStepClusters")
    ),
    FPix = cms.PSet(
        HitProducer = cms.string('siPixelRecHits'),
        TTRHBuilder = cms.string('WithTrackAngle'),
        skipClusters = cms.InputTag("pixelPairStepClusters")
    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(

    ),
    TIB = cms.PSet(

    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(

    ),
    layerList = cms.vstring(
        'BPix1+BPix2',
        'BPix1+BPix3',
        'BPix2+BPix3',
        'BPix1+FPix1_pos',
        'BPix1+FPix1_neg',
        'BPix2+FPix1_pos',
        'BPix2+FPix1_neg'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.pixelPairStepSeeds = cms.EDProducer("SeedCombiner",
    seedCollections = cms.VInputTag("pixelPairStepSeedsA", "pixelPairStepSeedsB")
)


process.pixelPairStepSeedsA = cms.EDProducer("SeedCreatorFromRegionConsecutiveHitsEDProducer",
    MinOneOverPtError = cms.double(1),
    OriginTransverseErrorMultiplier = cms.double(1),
    SeedComparitorPSet = cms.PSet(
        ClusterShapeCacheSrc = cms.InputTag("siPixelClusterShapeCache"),
        ClusterShapeHitFilterName = cms.string('ClusterShapeHitFilter'),
        ComponentName = cms.string('PixelClusterShapeSeedComparitor'),
        FilterAtHelixStage = cms.bool(True),
        FilterPixelHits = cms.bool(True),
        FilterStripHits = cms.bool(False)
    ),
    SeedMomentumForBOFF = cms.double(5),
    TTRHBuilder = cms.string('WithTrackAngle'),
    forceKinematicWithRegionDirection = cms.bool(False),
    magneticField = cms.string('ParabolicMf'),
    mightGet = cms.optional.untracked.vstring,
    propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    seedingHitSets = cms.InputTag("pixelPairStepHitDoublets")
)


process.pixelPairStepSeedsB = cms.EDProducer("SeedCreatorFromRegionConsecutiveHitsEDProducer",
    MinOneOverPtError = cms.double(1),
    OriginTransverseErrorMultiplier = cms.double(1),
    SeedComparitorPSet = cms.PSet(
        ClusterShapeCacheSrc = cms.InputTag("siPixelClusterShapeCache"),
        ClusterShapeHitFilterName = cms.string('ClusterShapeHitFilter'),
        ComponentName = cms.string('PixelClusterShapeSeedComparitor'),
        FilterAtHelixStage = cms.bool(True),
        FilterPixelHits = cms.bool(True),
        FilterStripHits = cms.bool(False)
    ),
    SeedMomentumForBOFF = cms.double(5),
    TTRHBuilder = cms.string('WithTrackAngle'),
    forceKinematicWithRegionDirection = cms.bool(False),
    magneticField = cms.string('ParabolicMf'),
    mightGet = cms.optional.untracked.vstring,
    propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    seedingHitSets = cms.InputTag("pixelPairStepHitDoubletsB")
)


process.pixelPairStepSelector = cms.EDProducer("MultiTrackSelector",
    GBRForestLabel = cms.string('MVASelectorIter2'),
    beamspot = cms.InputTag("offlineBeamSpot"),
    src = cms.InputTag("pixelPairStepTracks"),
    trackSelectors = cms.VPSet(
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(True),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(1.6),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(0.55, 4.0),
            d0_par2 = cms.vdouble(0.55, 4.0),
            dz_par1 = cms.vdouble(0.65, 4.0),
            dz_par2 = cms.vdouble(0.45, 4.0),
            keepAllTracks = cms.bool(False),
            maxNumberLostLayers = cms.uint32(999),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(0),
            minNumberLayers = cms.uint32(0),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(0),
            nSigmaZ = cms.double(4.0),
            name = cms.string('pixelPairStepLoose'),
            preFilterName = cms.string(''),
            qualityBit = cms.string('loose'),
            res_par = cms.vdouble(0.003, 0.01),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        ),
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(True),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(0.7),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(0.3, 4.0),
            d0_par2 = cms.vdouble(0.4, 4.0),
            dz_par1 = cms.vdouble(0.35, 4.0),
            dz_par2 = cms.vdouble(0.4, 4.0),
            keepAllTracks = cms.bool(True),
            maxNumberLostLayers = cms.uint32(2),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(3),
            minNumberLayers = cms.uint32(3),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(0),
            nSigmaZ = cms.double(4.0),
            name = cms.string('pixelPairStepTight'),
            preFilterName = cms.string('pixelPairStepLoose'),
            qualityBit = cms.string('tight'),
            res_par = cms.vdouble(0.003, 0.01),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        ),
        cms.PSet(
            applyAbsCutsIfNoPV = cms.bool(False),
            applyAdaptedPVCuts = cms.bool(True),
            chi2n_no1Dmod_par = cms.double(9999),
            chi2n_par = cms.double(0.7),
            copyExtras = cms.untracked.bool(True),
            copyTrajectories = cms.untracked.bool(False),
            d0_par1 = cms.vdouble(0.3, 4.0),
            d0_par2 = cms.vdouble(0.4, 4.0),
            dz_par1 = cms.vdouble(0.35, 4.0),
            dz_par2 = cms.vdouble(0.4, 4.0),
            keepAllTracks = cms.bool(True),
            maxNumberLostLayers = cms.uint32(2),
            max_d0 = cms.double(100.0),
            max_eta = cms.double(9999.0),
            max_lostHitFraction = cms.double(1.0),
            max_minMissHitOutOrIn = cms.int32(99),
            max_relpterr = cms.double(9999.0),
            max_z0 = cms.double(100.0),
            minHitsToBypassChecks = cms.uint32(20),
            minNumber3DLayers = cms.uint32(3),
            minNumberLayers = cms.uint32(3),
            min_eta = cms.double(-9999.0),
            min_nhits = cms.uint32(0),
            nSigmaZ = cms.double(4.0),
            name = cms.string('QualityMasks'),
            preFilterName = cms.string('pixelPairStepTight'),
            qualityBit = cms.string('highPurity'),
            res_par = cms.vdouble(0.003, 0.001),
            vertexCut = cms.string('ndof>=2&!isFake'),
            vtxNumber = cms.int32(-1)
        )
    ),
    useAnyMVA = cms.bool(True),
    useVertices = cms.bool(True),
    useVtxError = cms.bool(False),
    vertices = cms.InputTag("pixelVertices")
)


process.pixelPairStepTrackCandidates = cms.EDProducer("CkfTrackCandidateMaker",
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    RedundantSeedCleaner = cms.string('CachingSeedCleanerBySharedInput'),
    TrajectoryBuilderPSet = cms.PSet(
        refToPSet_ = cms.string('pixelPairStepTrajectoryBuilder')
    ),
    TrajectoryCleaner = cms.string('TrajectoryCleanerBySharedHits'),
    TransientInitialStateEstimatorParameters = cms.PSet(
        numberMeasurementsForFit = cms.int32(4),
        propagatorAlongTISE = cms.string('PropagatorWithMaterial'),
        propagatorOppositeTISE = cms.string('PropagatorWithMaterialOpposite')
    ),
    cleanTrajectoryAfterInOut = cms.bool(True),
    clustersToSkip = cms.InputTag("pixelPairStepClusters"),
    doSeedingRegionRebuilding = cms.bool(True),
    maxNSeeds = cms.uint32(500000),
    maxSeedsBeforeCleaning = cms.uint32(5000),
    mightGet = cms.optional.untracked.vstring,
    numHitsForSeedCleaner = cms.int32(50),
    onlyPixelHitsForSeedCleaner = cms.bool(True),
    phase2clustersToSkip = cms.InputTag(""),
    reverseTrajectories = cms.bool(False),
    src = cms.InputTag("pixelPairStepSeeds"),
    useHitsSplitting = cms.bool(True)
)


process.pixelPairStepTrackCandidatesMkFit = cms.EDProducer("MkFitProducer",
    backwardFitInCMSSW = cms.bool(False),
    buildingRoutine = cms.string('cloneEngine'),
    clustersToSkip = cms.InputTag("pixelPairStepClusters"),
    config = cms.ESInputTag("","pixelPairStepTrackCandidatesMkFitConfig"),
    eventOfHits = cms.InputTag("mkFitEventOfHits"),
    limitConcurrency = cms.untracked.bool(False),
    mightGet = cms.optional.untracked.vstring,
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutLoose')
    ),
    mkFitSilent = cms.untracked.bool(True),
    pixelHits = cms.InputTag("mkFitSiPixelHits"),
    removeDuplicates = cms.bool(True),
    seedCleaning = cms.bool(True),
    seeds = cms.InputTag("pixelPairStepTrackCandidatesMkFitSeeds"),
    stripHits = cms.InputTag("mkFitSiStripHits")
)


process.pixelPairStepTrackCandidatesMkFitSeeds = cms.EDProducer("MkFitSeedConverter",
    maxNSeeds = cms.uint32(500000),
    mightGet = cms.optional.untracked.vstring,
    seeds = cms.InputTag("pixelPairStepSeeds"),
    ttrhBuilder = cms.ESInputTag("","WithTrackAngle")
)


process.pixelPairStepTrackingRegions = cms.EDProducer("GlobalTrackingRegionWithVerticesEDProducer",
    RegionPSet = cms.PSet(
        VertexCollection = cms.InputTag("firstStepPrimaryVertices"),
        beamSpot = cms.InputTag("offlineBeamSpot"),
        fixedError = cms.double(0.03),
        halfLengthScaling4BigEvts = cms.bool(False),
        maxNVertices = cms.int32(5),
        maxPtMin = cms.double(1000),
        minHalfLength = cms.double(0),
        minOriginR = cms.double(0),
        nSigmaZ = cms.double(4),
        originRScaling4BigEvts = cms.bool(False),
        originRadius = cms.double(0.015),
        pixelClustersForScaling = cms.InputTag("siPixelClusters"),
        precise = cms.bool(True),
        ptMin = cms.double(0.6),
        ptMinScaling4BigEvts = cms.bool(False),
        scalingEndNPix = cms.double(1),
        scalingStartNPix = cms.double(0),
        sigmaZVertex = cms.double(3),
        useFakeVertices = cms.bool(False),
        useFixedError = cms.bool(True),
        useFoundVertices = cms.bool(True),
        useMultipleScattering = cms.bool(False)
    ),
    mightGet = cms.optional.untracked.vstring
)


process.pixelPairStepTrackingRegionsSeedLayersB = cms.EDProducer("PixelInactiveAreaTrackingRegionsSeedingLayersProducer",
    BPix = cms.PSet(
        HitProducer = cms.string('siPixelRecHits'),
        TTRHBuilder = cms.string('WithTrackAngle'),
        skipClusters = cms.InputTag("pixelPairStepClusters")
    ),
    FPix = cms.PSet(
        HitProducer = cms.string('siPixelRecHits'),
        TTRHBuilder = cms.string('WithTrackAngle'),
        skipClusters = cms.InputTag("pixelPairStepClusters")
    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    RegionPSet = cms.PSet(
        beamSpot = cms.InputTag("offlineBeamSpot"),
        deltaEta_Cand = cms.double(-1),
        deltaPhi_Cand = cms.double(-1),
        extraEta = cms.double(0),
        extraPhi = cms.double(0),
        input = cms.InputTag(""),
        maxNVertices = cms.int32(5),
        measurementTrackerName = cms.InputTag(""),
        nSigmaZBeamSpot = cms.double(4),
        nSigmaZVertex = cms.double(3),
        operationMode = cms.string('VerticesFixed'),
        originRadius = cms.double(0.015),
        precise = cms.bool(True),
        ptMin = cms.double(0.6),
        searchOpt = cms.bool(False),
        seedingMode = cms.string('Global'),
        vertexCollection = cms.InputTag("firstStepPrimaryVertices"),
        whereToUseMeasurementTracker = cms.string('Never'),
        zErrorBeamSpot = cms.double(24.2),
        zErrorVertex = cms.double(0.03)
    ),
    TEC = cms.PSet(

    ),
    TIB = cms.PSet(

    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(

    ),
    badPixelFEDChannelCollectionLabels = cms.VInputTag("siPixelDigis"),
    createPlottingFiles = cms.untracked.bool(False),
    debug = cms.untracked.bool(False),
    ignoreSingleFPixPanelModules = cms.bool(True),
    inactivePixelDetectorLabels = cms.VInputTag("siPixelDigis"),
    layerList = cms.vstring(
        'BPix1+BPix4',
        'BPix2+BPix4',
        'BPix3+BPix4',
        'BPix1+FPix2_pos',
        'BPix1+FPix2_neg',
        'BPix1+FPix3_pos',
        'BPix1+FPix3_neg',
        'BPix2+FPix2_pos',
        'BPix2+FPix2_neg',
        'BPix3+FPix1_pos',
        'BPix3+FPix1_neg',
        'FPix1_pos+FPix2_pos',
        'FPix1_neg+FPix2_neg',
        'FPix1_pos+FPix3_pos',
        'FPix1_neg+FPix3_neg',
        'FPix2_pos+FPix3_pos',
        'FPix2_neg+FPix3_neg'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.pixelPairStepTracks = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('pixelPairStep'),
    Fitter = cms.string('FlexibleKFFittingSmoother'),
    GeometricInnerState = cms.bool(False),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    Propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    TTRHBuilder = cms.string('WithAngleAndTemplateWithoutProbQ'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("pixelPairStepTrackCandidates"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(True)
)


process.preDuplicateMergingGeneralTracks = cms.EDProducer("TrackCollectionMerger",
    allowFirstHitShare = cms.bool(True),
    copyExtras = cms.untracked.bool(True),
    copyTrajectories = cms.untracked.bool(False),
    enableMerging = cms.bool(True),
    foundHitBonus = cms.double(100.0),
    inputClassifiers = cms.vstring(
        'earlyGeneralTracks',
        'muonSeededTracksInOutClassifier',
        'muonSeededTracksOutInClassifier'
    ),
    lostHitPenalty = cms.double(1.0),
    mightGet = cms.optional.untracked.vstring,
    minQuality = cms.string('loose'),
    minShareHits = cms.uint32(2),
    shareFrac = cms.double(0.19),
    trackAlgoPriorityOrder = cms.string('trackAlgoPriorityOrder'),
    trackProducers = cms.VInputTag("earlyGeneralTracks", "muonSeededTracksInOut", "muonSeededTracksOutIn")
)


process.recoTauAK4PFJets08Region = cms.EDProducer("RecoTauJetRegionProducer",
    deltaR = cms.double(0.8),
    maxJetAbsEta = cms.double(2.7),
    minJetPt = cms.double(14.0),
    pfCandAssocMapSrc = cms.InputTag(""),
    pfCandSrc = cms.InputTag("particleFlow"),
    src = cms.InputTag("ak4PFJets"),
    verbosity = cms.int32(0)
)


process.recoTauDiscriminantCutMultiplexerDefault = cms.EDProducer("RecoTauDiscriminantCutMultiplexer",
    PFTauProducer = cms.InputTag("fixme"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('AND'),
        decayMode = cms.PSet(
            Producer = cms.InputTag("fixme"),
            cut = cms.double(0)
        ),
        leadTrack = cms.PSet(
            Producer = cms.InputTag("fixme"),
            cut = cms.double(0)
        )
    ),
    inputFileName = cms.FileInPath('RecoTauTag/RecoTau/data/emptyMVAinputFile'),
    loadMVAfromDB = cms.bool(True),
    mapping = cms.VPSet(cms.PSet(
        category = cms.uint32(0),
        cut = cms.string('fixme')
    )),
    mightGet = cms.optional.untracked.vstring,
    mvaOutput_normalization = cms.string(''),
    rawValues = cms.vstring('discriminator'),
    toMultiplex = cms.InputTag("fixme"),
    verbosity = cms.int32(0),
    workingPoints = cms.vdouble(0)
)


process.regionalCosmicCkfTrackCandidates = cms.EDProducer("CkfTrackCandidateMaker",
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('CosmicNavigationSchool'),
    RedundantSeedCleaner = cms.string('CachingSeedCleanerBySharedInput'),
    TrajectoryBuilderPSet = cms.PSet(
        refToPSet_ = cms.string('GroupedCkfTrajectoryBuilderP5')
    ),
    TrajectoryCleaner = cms.string('TrajectoryCleanerBySharedHits'),
    TransientInitialStateEstimatorParameters = cms.PSet(
        numberMeasurementsForFit = cms.int32(4),
        propagatorAlongTISE = cms.string('PropagatorWithMaterial'),
        propagatorOppositeTISE = cms.string('PropagatorWithMaterialOpposite')
    ),
    cleanTrajectoryAfterInOut = cms.bool(True),
    clustersToSkip = cms.InputTag(""),
    doSeedingRegionRebuilding = cms.bool(True),
    maxNSeeds = cms.uint32(500000),
    maxSeedsBeforeCleaning = cms.uint32(5000),
    mightGet = cms.optional.untracked.vstring,
    numHitsForSeedCleaner = cms.int32(4),
    onlyPixelHitsForSeedCleaner = cms.bool(False),
    phase2clustersToSkip = cms.InputTag(""),
    reverseTrajectories = cms.bool(False),
    src = cms.InputTag("regionalCosmicTrackerSeeds"),
    useHitsSplitting = cms.bool(True)
)


process.regionalCosmicTrackerSeedingLayers = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(

    ),
    FPix = cms.PSet(

    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        maxRing = cms.int32(7),
        minRing = cms.int32(6),
        useRingSlector = cms.bool(False)
    ),
    TIB = cms.PSet(

    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        )
    ),
    layerList = cms.vstring(
        'TOB6+TOB5',
        'TOB6+TOB4',
        'TOB6+TOB3',
        'TOB5+TOB4',
        'TOB5+TOB3',
        'TOB4+TOB3',
        'TEC1_neg+TOB6',
        'TEC1_neg+TOB5',
        'TEC1_neg+TOB4',
        'TEC1_pos+TOB6',
        'TEC1_pos+TOB5',
        'TEC1_pos+TOB4'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.regionalCosmicTrackerSeeds = cms.EDProducer("SeedGeneratorFromRegionHitsEDProducer",
    ClusterCheckPSet = cms.PSet(
        ClusterCollectionLabel = cms.InputTag("siStripClusters"),
        MaxNumberOfCosmicClusters = cms.uint32(10000),
        MaxNumberOfPixelClusters = cms.uint32(10000),
        PixelClusterCollectionLabel = cms.InputTag("siPixelClusters"),
        doClusterCheck = cms.bool(False)
    ),
    OrderedHitsFactoryPSet = cms.PSet(
        ComponentName = cms.string('GenericPairGenerator'),
        LayerSrc = cms.InputTag("regionalCosmicTrackerSeedingLayers")
    ),
    RegionFactoryPSet = cms.PSet(
        CollectionsPSet = cms.PSet(
            recoL2MuonsCollection = cms.InputTag(""),
            recoMuonsCollection = cms.InputTag(""),
            recoTrackMuonsCollection = cms.InputTag("cosmicMuons")
        ),
        ComponentName = cms.string('CosmicRegionalSeedGenerator'),
        RegionInJetsCheckPSet = cms.PSet(
            deltaRExclusionSize = cms.double(0.3),
            doJetsExclusionCheck = cms.bool(True),
            jetsPtMin = cms.double(5),
            recoCaloJetsCollection = cms.InputTag("ak4CaloJets")
        ),
        RegionPSet = cms.PSet(
            deltaEtaRegion = cms.double(0.1),
            deltaPhiRegion = cms.double(0.1),
            measurementTrackerName = cms.string(''),
            precise = cms.bool(True),
            ptMin = cms.double(1.0),
            rVertex = cms.double(5),
            zVertex = cms.double(5)
        ),
        ToolsPSet = cms.PSet(
            regionBase = cms.string('seedOnCosmicMuon'),
            thePropagatorName = cms.string('AnalyticalPropagator')
        )
    ),
    RegionInJetsCheckPSet = cms.PSet(
        doJetsExclusionCheck = cms.bool(False)
    ),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('none')
    ),
    SeedCreatorPSet = cms.PSet(
        ComponentName = cms.string('CosmicSeedCreator'),
        MinOneOverPtError = cms.double(1.0),
        OriginTransverseErrorMultiplier = cms.double(1.0),
        SeedMomentumForBOFF = cms.double(5.0),
        TTRHBuilder = cms.string('WithTrackAngle'),
        forceKinematicWithRegionDirection = cms.bool(False),
        magneticField = cms.string('ParabolicMf'),
        maxseeds = cms.int32(10000),
        propagator = cms.string('PropagatorWithMaterialParabolicMf')
    )
)


process.regionalCosmicTracks = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('ctf'),
    Fitter = cms.string('FittingSmootherRKP5'),
    GeometricInnerState = cms.bool(True),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    Propagator = cms.string('RungeKuttaTrackerPropagator'),
    SimpleMagneticField = cms.string(''),
    TTRHBuilder = cms.string('WithAngleAndTemplate'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("regionalCosmicCkfTrackCandidates"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(False)
)


process.seedClusterRemover = cms.EDProducer("SeedClusterRemover",
    Common = cms.PSet(
        maxChi2 = cms.double(9.0)
    ),
    pixelClusters = cms.InputTag("siPixelClusters"),
    stripClusters = cms.InputTag("siStripClusters"),
    trajectories = cms.InputTag("initialStepSeeds")
)


process.seedClusterRemoverPhase2 = cms.EDProducer("SeedClusterRemoverPhase2",
    phase2OTClusters = cms.InputTag("siPhase2Clusters"),
    pixelClusters = cms.InputTag("siPixelClusters"),
    trajectories = cms.InputTag("initialStepSeeds")
)


process.seedGeneratorFromRegionHitsEDProducer = cms.EDProducer("SeedGeneratorFromRegionHitsEDProducer",
    ClusterCheckPSet = cms.PSet(
        ClusterCollectionLabel = cms.InputTag("siStripClusters"),
        MaxNumberOfCosmicClusters = cms.uint32(400000),
        MaxNumberOfPixelClusters = cms.uint32(40000),
        PixelClusterCollectionLabel = cms.InputTag("siPixelClusters"),
        cut = cms.string('strip < 400000 && pixel < 40000 && (strip < 50000 + 10*pixel) && (pixel < 5000 + 0.1*strip)'),
        doClusterCheck = cms.bool(False)
    ),
    OrderedHitsFactoryPSet = cms.PSet(
        ComponentName = cms.string(''),
        SeedingLayers = cms.InputTag(""),
        maxElement = cms.uint32(1000000)
    ),
    RegionFactoryPSet = cms.PSet(
        ComponentName = cms.string('GlobalRegionProducerFromBeamSpot'),
        RegionPSet = cms.PSet(
            beamSpot = cms.InputTag("offlineBeamSpot"),
            originHalfLength = cms.double(21.2),
            originRadius = cms.double(0.2),
            precise = cms.bool(True),
            ptMin = cms.double(0.9),
            useMultipleScattering = cms.bool(False)
        )
    ),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('none')
    ),
    SeedCreatorPSet = cms.PSet(
        ComponentName = cms.string('SeedFromConsecutiveHitsCreator'),
        MinOneOverPtError = cms.double(1.0),
        OriginTransverseErrorMultiplier = cms.double(1.0),
        SeedMomentumForBOFF = cms.double(5.0),
        TTRHBuilder = cms.string('WithTrackAngle'),
        forceKinematicWithRegionDirection = cms.bool(False),
        magneticField = cms.string('ParabolicMf'),
        propagator = cms.string('PropagatorWithMaterialParabolicMf')
    )
)


process.seedingLayersEDProducer = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(

    ),
    FPix = cms.PSet(

    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(

    ),
    TIB = cms.PSet(

    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(

    ),
    layerList = cms.vstring(),
    mightGet = cms.optional.untracked.vstring
)


process.siPhase2Clusters = cms.EDProducer("Phase2TrackerClusterizer",
    maxClusterSize = cms.uint32(0),
    maxNumberClusters = cms.uint32(0),
    src = cms.InputTag("mix","Tracker")
)


process.siPhase2VectorHits = cms.EDProducer("VectorHitBuilderEDProducer",
    Algorithm = cms.ESInputTag("","SiPhase2VectorHitMatcher"),
    BarrelCut = cms.vdouble(
        0, 0.05, 0.06, 0.08, 0.09,
        0.12, 0.2
    ),
    CPE = cms.ESInputTag("phase2StripCPEESProducer","Phase2StripCPE"),
    Clusters = cms.string('siPhase2Clusters'),
    EndcapCut = cms.vdouble(
        0, 0.1, 0.1, 0.1, 0.1,
        0.1
    ),
    maxVectorHits = cms.int32(999999999),
    maxVectorHitsInAStack = cms.int32(999),
    mightGet = cms.optional.untracked.vstring,
    offlinestubs = cms.string('vectorHits')
)


process.siPixelClusterShapeCache = cms.EDProducer("SiPixelClusterShapeCacheProducer",
    mightGet = cms.optional.untracked.vstring,
    onDemand = cms.bool(False),
    src = cms.InputTag("siPixelClusters")
)


process.siPixelClusters = cms.EDProducer("JetCoreClusterSplitter",
    centralMIPCharge = cms.double(26000),
    chargeFractionMin = cms.double(2.0),
    chargePerUnit = cms.double(2000),
    cores = cms.InputTag("jetsForCoreTrackingPreSplitting"),
    deltaRmax = cms.double(0.05),
    forceXError = cms.double(100),
    forceYError = cms.double(150),
    fractionalWidth = cms.double(0.4),
    pixelCPE = cms.string('PixelCPEGeneric'),
    pixelClusters = cms.InputTag("siPixelClustersPreSplitting"),
    ptMin = cms.double(200),
    verbose = cms.bool(False),
    vertices = cms.InputTag("firstStepPrimaryVerticesPreSplitting")
)


process.siPixelClustersBottom = cms.EDProducer("PixelClusterSelectorTopBottom",
    label = cms.InputTag("siPixelClusters"),
    y = cms.double(-1)
)


process.siPixelClustersPreSplittingCUDA = cms.EDProducer("SiPixelRawToClusterCUDA",
    CablingMapLabel = cms.string(''),
    IncludeErrors = cms.bool(True),
    InputLabel = cms.InputTag("rawDataCollector"),
    Regions = cms.PSet(
        beamSpot = cms.optional.InputTag,
        deltaPhi = cms.optional.vdouble,
        inputs = cms.optional.VInputTag,
        maxZ = cms.optional.vdouble
    ),
    UseQualityInfo = cms.bool(False),
    clusterThreshold_layer1 = cms.int32(4000),
    clusterThreshold_otherLayers = cms.int32(4000),
    isRun2 = cms.bool(False),
    mightGet = cms.optional.untracked.vstring
)


process.siPixelClustersTop = cms.EDProducer("PixelClusterSelectorTopBottom",
    label = cms.InputTag("siPixelClusters"),
    y = cms.double(1)
)


process.siPixelDigisClustersPreSplitting = cms.EDProducer("SiPixelDigisClustersFromSoAPhase1",
    clusterThreshold_layer1 = cms.int32(4000),
    clusterThreshold_otherLayers = cms.int32(4000),
    mightGet = cms.optional.untracked.vstring,
    produceDigis = cms.bool(True),
    src = cms.InputTag("siPixelDigisSoA"),
    storeDigis = cms.bool(True)
)


process.siPixelDigisMorphed = cms.EDProducer("SiPixelDigiMorphing",
    fakeAdc = cms.uint32(100),
    iters = cms.int32(1),
    kernel1 = cms.vint32(7, 7, 7),
    kernel2 = cms.vint32(2, 7, 2),
    mightGet = cms.optional.untracked.vstring,
    ncols = cms.int32(416),
    nrocs = cms.int32(8),
    nrows = cms.int32(160),
    src = cms.InputTag("siPixelDigis")
)


process.siPixelDigisPhase2SoA = cms.EDProducer("SiPixelDigisSoAFromCUDA",
    mightGet = cms.optional.untracked.vstring,
    src = cms.InputTag("siPixelClustersPreSplittingCUDA")
)


process.siPixelRecHits = cms.EDProducer("SiPixelRecHitConverter",
    CPE = cms.string('PixelCPEGeneric'),
    VerboseLevel = cms.untracked.int32(0),
    src = cms.InputTag("siPixelClusters")
)


process.siPixelRecHitsBottom = cms.EDProducer("SiPixelRecHitConverter",
    CPE = cms.string('PixelCPEGeneric'),
    VerboseLevel = cms.untracked.int32(0),
    src = cms.InputTag("siPixelClustersBottom")
)


process.siPixelRecHitsPreSplittingCPU = cms.EDProducer("SiPixelRecHitSoAFromLegacyPhase1",
    CPE = cms.string('PixelCPEFast'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    convertToLegacy = cms.bool(True),
    mightGet = cms.optional.untracked.vstring,
    src = cms.InputTag("siPixelClustersPreSplitting")
)


process.siPixelRecHitsPreSplittingCUDA = cms.EDProducer("SiPixelRecHitCUDAPhase1",
    CPE = cms.string('PixelCPEFast'),
    beamSpot = cms.InputTag("offlineBeamSpotToCUDA"),
    mightGet = cms.optional.untracked.vstring,
    src = cms.InputTag("siPixelClustersPreSplittingCUDA")
)


process.siPixelRecHitsTop = cms.EDProducer("SiPixelRecHitConverter",
    CPE = cms.string('PixelCPEGeneric'),
    VerboseLevel = cms.untracked.int32(0),
    src = cms.InputTag("siPixelClustersTop")
)


process.siStripClusters = cms.EDProducer("SiStripClusterizer",
    Clusterizer = cms.PSet(
        Algorithm = cms.string('ThreeThresholdAlgorithm'),
        ChannelThreshold = cms.double(2.0),
        ClusterThreshold = cms.double(5.0),
        ConditionsLabel = cms.string(''),
        MaxAdjacentBad = cms.uint32(0),
        MaxSequentialBad = cms.uint32(1),
        MaxSequentialHoles = cms.uint32(0),
        RemoveApvShots = cms.bool(True),
        SeedThreshold = cms.double(3.0),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        )
    ),
    DigiProducersList = cms.VInputTag(cms.InputTag("siStripDigis","ZeroSuppressed"), cms.InputTag("siStripZeroSuppression","VirginRaw"), cms.InputTag("siStripZeroSuppression","ProcessedRaw"), cms.InputTag("siStripZeroSuppression","ScopeMode"))
)


process.siStripClustersBottom = cms.EDProducer("StripClusterSelectorTopBottom",
    label = cms.InputTag("siStripClusters"),
    y = cms.double(-1)
)


process.siStripClustersTop = cms.EDProducer("StripClusterSelectorTopBottom",
    label = cms.InputTag("siStripClusters"),
    y = cms.double(1)
)


process.siStripMatchedRecHits = cms.EDProducer("SiStripRecHitConverter",
    ClusterProducer = cms.InputTag("siStripClusters"),
    MaskBadAPVFibers = cms.bool(False),
    Matcher = cms.ESInputTag("SiStripRecHitMatcherESProducer","StandardMatcher"),
    StripCPE = cms.ESInputTag("StripCPEfromTrackAngleESProducer","StripCPEfromTrackAngle"),
    VerbosityLevel = cms.optional.untracked.int32,
    doMatching = cms.bool(True),
    matchedRecHits = cms.string('matchedRecHit'),
    mightGet = cms.optional.untracked.vstring,
    rphiRecHits = cms.string('rphiRecHit'),
    siStripQualityLabel = cms.ESInputTag("",""),
    stereoRecHits = cms.string('stereoRecHit'),
    useSiStripQuality = cms.bool(False)
)


process.siStripMatchedRecHitsBottom = cms.EDProducer("SiStripRecHitConverter",
    ClusterProducer = cms.InputTag("siStripClustersBottom"),
    MaskBadAPVFibers = cms.bool(False),
    Matcher = cms.ESInputTag("SiStripRecHitMatcherESProducer","StandardMatcher"),
    StripCPE = cms.ESInputTag("StripCPEfromTrackAngleESProducer","StripCPEfromTrackAngle"),
    VerbosityLevel = cms.optional.untracked.int32,
    doMatching = cms.bool(True),
    matchedRecHits = cms.string('matchedRecHit'),
    mightGet = cms.optional.untracked.vstring,
    rphiRecHits = cms.string('rphiRecHit'),
    siStripQualityLabel = cms.ESInputTag("",""),
    stereoRecHits = cms.string('stereoRecHit'),
    useSiStripQuality = cms.bool(False)
)


process.siStripMatchedRecHitsTop = cms.EDProducer("SiStripRecHitConverter",
    ClusterProducer = cms.InputTag("siStripClustersTop"),
    MaskBadAPVFibers = cms.bool(False),
    Matcher = cms.ESInputTag("SiStripRecHitMatcherESProducer","StandardMatcher"),
    StripCPE = cms.ESInputTag("StripCPEfromTrackAngleESProducer","StripCPEfromTrackAngle"),
    VerbosityLevel = cms.optional.untracked.int32,
    doMatching = cms.bool(True),
    matchedRecHits = cms.string('matchedRecHit'),
    mightGet = cms.optional.untracked.vstring,
    rphiRecHits = cms.string('rphiRecHit'),
    siStripQualityLabel = cms.ESInputTag("",""),
    stereoRecHits = cms.string('stereoRecHit'),
    useSiStripQuality = cms.bool(False)
)


process.siStripZeroSuppression = cms.EDProducer("SiStripZeroSuppression",
    Algorithms = cms.PSet(
        APVInspectMode = cms.string('BaselineFollower'),
        APVRestoreMode = cms.string('BaselineFollower'),
        ApplyBaselineCleaner = cms.bool(True),
        ApplyBaselineRejection = cms.bool(True),
        CleaningSequence = cms.uint32(1),
        CommonModeNoiseSubtractionMode = cms.string('IteratedMedian'),
        CutToAvoidSignal = cms.double(2.0),
        DeltaCMThreshold = cms.uint32(20),
        Deviation = cms.uint32(25),
        ForceNoRestore = cms.bool(False),
        Fraction = cms.double(0.2),
        Iterations = cms.int32(3),
        MeanCM = cms.int32(0),
        PedestalSubtractionFedMode = cms.bool(False),
        SiStripFedZeroSuppressionMode = cms.uint32(4),
        TruncateInSuppressor = cms.bool(True),
        Use10bitsTruncation = cms.bool(False),
        consecThreshold = cms.uint32(5),
        discontinuityThreshold = cms.int32(12),
        distortionThreshold = cms.uint32(20),
        doAPVRestore = cms.bool(True),
        filteredBaselineDerivativeSumSquare = cms.double(30),
        filteredBaselineMax = cms.double(6),
        hitStripThreshold = cms.uint32(40),
        lastGradient = cms.int32(10),
        minStripsToFit = cms.uint32(4),
        nSaturatedStrip = cms.uint32(2),
        nSigmaNoiseDerTh = cms.uint32(4),
        nSmooth = cms.uint32(9),
        restoreThreshold = cms.double(0.5),
        sizeWindow = cms.int32(1),
        slopeX = cms.int32(3),
        slopeY = cms.int32(4),
        useCMMeanMap = cms.bool(False),
        useRealMeanCM = cms.bool(False),
        widthCluster = cms.int32(64)
    ),
    RawDigiProducersList = cms.VInputTag(cms.InputTag("siStripDigis","VirginRaw"), cms.InputTag("siStripDigis","ProcessedRaw"), cms.InputTag("siStripDigis","ScopeMode")),
    fixCM = cms.bool(False),
    produceBaselinePoints = cms.bool(False),
    produceCalculatedBaseline = cms.bool(False),
    produceHybridFormat = cms.bool(False),
    produceRawDigis = cms.bool(True),
    storeCM = cms.bool(True),
    storeInZScollBadAPV = cms.bool(True)
)


process.siStripZeroSuppressionHLT = cms.EDProducer("SiStripZeroSuppression",
    Algorithms = cms.PSet(
        APVInspectMode = cms.string('BaselineFollower'),
        APVRestoreMode = cms.string('BaselineFollower'),
        ApplyBaselineCleaner = cms.bool(True),
        ApplyBaselineRejection = cms.bool(True),
        CleaningSequence = cms.uint32(1),
        CommonModeNoiseSubtractionMode = cms.string('IteratedMedian'),
        CutToAvoidSignal = cms.double(2.0),
        DeltaCMThreshold = cms.uint32(20),
        Deviation = cms.uint32(25),
        ForceNoRestore = cms.bool(False),
        Fraction = cms.double(0.2),
        Iterations = cms.int32(3),
        MeanCM = cms.int32(0),
        PedestalSubtractionFedMode = cms.bool(False),
        SiStripFedZeroSuppressionMode = cms.uint32(4),
        TruncateInSuppressor = cms.bool(True),
        Use10bitsTruncation = cms.bool(False),
        consecThreshold = cms.uint32(5),
        discontinuityThreshold = cms.int32(12),
        distortionThreshold = cms.uint32(20),
        doAPVRestore = cms.bool(True),
        filteredBaselineDerivativeSumSquare = cms.double(30),
        filteredBaselineMax = cms.double(6),
        hitStripThreshold = cms.uint32(40),
        lastGradient = cms.int32(10),
        minStripsToFit = cms.uint32(4),
        nSaturatedStrip = cms.uint32(2),
        nSigmaNoiseDerTh = cms.uint32(4),
        nSmooth = cms.uint32(9),
        restoreThreshold = cms.double(0.5),
        sizeWindow = cms.int32(1),
        slopeX = cms.int32(3),
        slopeY = cms.int32(4),
        useCMMeanMap = cms.bool(False),
        useRealMeanCM = cms.bool(False),
        widthCluster = cms.int32(64)
    ),
    RawDigiProducersList = cms.VInputTag("siStripDigisHLT:VirginRaw", "siStripDigisHLT:ProcessedRaw", "siStripDigisHLT:ScopeMode"),
    fixCM = cms.bool(False),
    produceBaselinePoints = cms.bool(False),
    produceCalculatedBaseline = cms.bool(False),
    produceHybridFormat = cms.bool(False),
    produceRawDigis = cms.bool(True),
    storeCM = cms.bool(True),
    storeInZScollBadAPV = cms.bool(True)
)


process.simpleCosmicBONSeedingLayers = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(

    ),
    FPix = cms.PSet(

    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit")
    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit")
    ),
    TEC = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        maxRing = cms.int32(7),
        minRing = cms.int32(5),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
        useRingSlector = cms.bool(False),
        useSimpleRphiHitsCleaner = cms.bool(False)
    ),
    TIB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit")
    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit")
    ),
    layerList = cms.vstring(
        'MTOB4+MTOB5+MTOB6',
        'MTOB3+MTOB5+MTOB6',
        'MTOB3+MTOB4+MTOB5',
        'MTOB3+MTOB4+MTOB6',
        'TOB2+MTOB4+MTOB5',
        'TOB2+MTOB3+MTOB5',
        'TEC7_pos+TEC8_pos+TEC9_pos',
        'TEC6_pos+TEC7_pos+TEC8_pos',
        'TEC5_pos+TEC6_pos+TEC7_pos',
        'TEC4_pos+TEC5_pos+TEC6_pos',
        'TEC3_pos+TEC4_pos+TEC5_pos',
        'TEC2_pos+TEC3_pos+TEC4_pos',
        'TEC1_pos+TEC2_pos+TEC3_pos',
        'TEC7_neg+TEC8_neg+TEC9_neg',
        'TEC6_neg+TEC7_neg+TEC8_neg',
        'TEC5_neg+TEC6_neg+TEC7_neg',
        'TEC4_neg+TEC5_neg+TEC6_neg',
        'TEC3_neg+TEC4_neg+TEC5_neg',
        'TEC2_neg+TEC3_neg+TEC4_neg',
        'TEC1_neg+TEC2_neg+TEC3_neg',
        'TEC6_pos+TEC8_pos+TEC9_pos',
        'TEC5_pos+TEC7_pos+TEC8_pos',
        'TEC4_pos+TEC6_pos+TEC7_pos',
        'TEC3_pos+TEC5_pos+TEC6_pos',
        'TEC2_pos+TEC4_pos+TEC5_pos',
        'TEC1_pos+TEC3_pos+TEC4_pos',
        'TEC6_pos+TEC7_pos+TEC9_pos',
        'TEC5_pos+TEC6_pos+TEC8_pos',
        'TEC4_pos+TEC5_pos+TEC7_pos',
        'TEC3_pos+TEC4_pos+TEC6_pos',
        'TEC2_pos+TEC3_pos+TEC5_pos',
        'TEC1_pos+TEC2_pos+TEC4_pos',
        'TEC6_neg+TEC8_neg+TEC9_neg',
        'TEC5_neg+TEC7_neg+TEC8_neg',
        'TEC4_neg+TEC6_neg+TEC7_neg',
        'TEC3_neg+TEC5_neg+TEC6_neg',
        'TEC2_neg+TEC4_neg+TEC5_neg',
        'TEC1_neg+TEC3_neg+TEC4_neg',
        'TEC6_neg+TEC7_neg+TEC9_neg',
        'TEC5_neg+TEC6_neg+TEC8_neg',
        'TEC4_neg+TEC5_neg+TEC7_neg',
        'TEC3_neg+TEC4_neg+TEC6_neg',
        'TEC2_neg+TEC3_neg+TEC5_neg',
        'TEC1_neg+TEC2_neg+TEC4_neg',
        'MTOB6+TEC1_pos+TEC2_pos',
        'MTOB6+TEC1_neg+TEC2_neg',
        'MTOB6+MTOB5+TEC1_pos',
        'MTOB6+MTOB5+TEC1_neg'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.simpleCosmicBONSeedingLayersBottom = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(

    ),
    FPix = cms.PSet(

    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        rphiRecHits = cms.InputTag("siStripMatchedRecHitsBottom","rphiRecHit")
    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        rphiRecHits = cms.InputTag("siStripMatchedRecHitsBottom","rphiRecHit")
    ),
    TEC = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHitsBottom","matchedRecHit"),
        maxRing = cms.int32(7),
        minRing = cms.int32(5),
        rphiRecHits = cms.InputTag("siStripMatchedRecHitsBottom","rphiRecHit"),
        useRingSlector = cms.bool(False),
        useSimpleRphiHitsCleaner = cms.bool(False)
    ),
    TIB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHitsBottom","matchedRecHit")
    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHitsBottom","matchedRecHit")
    ),
    layerList = cms.vstring(
        'MTOB4+MTOB5+MTOB6',
        'MTOB3+MTOB5+MTOB6',
        'MTOB3+MTOB4+MTOB5',
        'MTOB3+MTOB4+MTOB6',
        'TOB2+MTOB4+MTOB5',
        'TOB2+MTOB3+MTOB5',
        'TEC7_pos+TEC8_pos+TEC9_pos',
        'TEC6_pos+TEC7_pos+TEC8_pos',
        'TEC5_pos+TEC6_pos+TEC7_pos',
        'TEC4_pos+TEC5_pos+TEC6_pos',
        'TEC3_pos+TEC4_pos+TEC5_pos',
        'TEC2_pos+TEC3_pos+TEC4_pos',
        'TEC1_pos+TEC2_pos+TEC3_pos',
        'TEC7_neg+TEC8_neg+TEC9_neg',
        'TEC6_neg+TEC7_neg+TEC8_neg',
        'TEC5_neg+TEC6_neg+TEC7_neg',
        'TEC4_neg+TEC5_neg+TEC6_neg',
        'TEC3_neg+TEC4_neg+TEC5_neg',
        'TEC2_neg+TEC3_neg+TEC4_neg',
        'TEC1_neg+TEC2_neg+TEC3_neg',
        'TEC6_pos+TEC8_pos+TEC9_pos',
        'TEC5_pos+TEC7_pos+TEC8_pos',
        'TEC4_pos+TEC6_pos+TEC7_pos',
        'TEC3_pos+TEC5_pos+TEC6_pos',
        'TEC2_pos+TEC4_pos+TEC5_pos',
        'TEC1_pos+TEC3_pos+TEC4_pos',
        'TEC6_pos+TEC7_pos+TEC9_pos',
        'TEC5_pos+TEC6_pos+TEC8_pos',
        'TEC4_pos+TEC5_pos+TEC7_pos',
        'TEC3_pos+TEC4_pos+TEC6_pos',
        'TEC2_pos+TEC3_pos+TEC5_pos',
        'TEC1_pos+TEC2_pos+TEC4_pos',
        'TEC6_neg+TEC8_neg+TEC9_neg',
        'TEC5_neg+TEC7_neg+TEC8_neg',
        'TEC4_neg+TEC6_neg+TEC7_neg',
        'TEC3_neg+TEC5_neg+TEC6_neg',
        'TEC2_neg+TEC4_neg+TEC5_neg',
        'TEC1_neg+TEC3_neg+TEC4_neg',
        'TEC6_neg+TEC7_neg+TEC9_neg',
        'TEC5_neg+TEC6_neg+TEC8_neg',
        'TEC4_neg+TEC5_neg+TEC7_neg',
        'TEC3_neg+TEC4_neg+TEC6_neg',
        'TEC2_neg+TEC3_neg+TEC5_neg',
        'TEC1_neg+TEC2_neg+TEC4_neg',
        'MTOB6+TEC1_pos+TEC2_pos',
        'MTOB6+TEC1_neg+TEC2_neg',
        'MTOB6+MTOB5+TEC1_pos',
        'MTOB6+MTOB5+TEC1_neg'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.simpleCosmicBONSeedingLayersTop = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(

    ),
    FPix = cms.PSet(

    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        rphiRecHits = cms.InputTag("siStripMatchedRecHitsTop","rphiRecHit")
    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        rphiRecHits = cms.InputTag("siStripMatchedRecHitsTop","rphiRecHit")
    ),
    TEC = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHitsTop","matchedRecHit"),
        maxRing = cms.int32(7),
        minRing = cms.int32(5),
        rphiRecHits = cms.InputTag("siStripMatchedRecHitsTop","rphiRecHit"),
        useRingSlector = cms.bool(False),
        useSimpleRphiHitsCleaner = cms.bool(False)
    ),
    TIB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHitsTop","matchedRecHit")
    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHitsTop","matchedRecHit")
    ),
    layerList = cms.vstring(
        'MTOB4+MTOB5+MTOB6',
        'MTOB3+MTOB5+MTOB6',
        'MTOB3+MTOB4+MTOB5',
        'MTOB3+MTOB4+MTOB6',
        'TOB2+MTOB4+MTOB5',
        'TOB2+MTOB3+MTOB5',
        'TEC7_pos+TEC8_pos+TEC9_pos',
        'TEC6_pos+TEC7_pos+TEC8_pos',
        'TEC5_pos+TEC6_pos+TEC7_pos',
        'TEC4_pos+TEC5_pos+TEC6_pos',
        'TEC3_pos+TEC4_pos+TEC5_pos',
        'TEC2_pos+TEC3_pos+TEC4_pos',
        'TEC1_pos+TEC2_pos+TEC3_pos',
        'TEC7_neg+TEC8_neg+TEC9_neg',
        'TEC6_neg+TEC7_neg+TEC8_neg',
        'TEC5_neg+TEC6_neg+TEC7_neg',
        'TEC4_neg+TEC5_neg+TEC6_neg',
        'TEC3_neg+TEC4_neg+TEC5_neg',
        'TEC2_neg+TEC3_neg+TEC4_neg',
        'TEC1_neg+TEC2_neg+TEC3_neg',
        'TEC6_pos+TEC8_pos+TEC9_pos',
        'TEC5_pos+TEC7_pos+TEC8_pos',
        'TEC4_pos+TEC6_pos+TEC7_pos',
        'TEC3_pos+TEC5_pos+TEC6_pos',
        'TEC2_pos+TEC4_pos+TEC5_pos',
        'TEC1_pos+TEC3_pos+TEC4_pos',
        'TEC6_pos+TEC7_pos+TEC9_pos',
        'TEC5_pos+TEC6_pos+TEC8_pos',
        'TEC4_pos+TEC5_pos+TEC7_pos',
        'TEC3_pos+TEC4_pos+TEC6_pos',
        'TEC2_pos+TEC3_pos+TEC5_pos',
        'TEC1_pos+TEC2_pos+TEC4_pos',
        'TEC6_neg+TEC8_neg+TEC9_neg',
        'TEC5_neg+TEC7_neg+TEC8_neg',
        'TEC4_neg+TEC6_neg+TEC7_neg',
        'TEC3_neg+TEC5_neg+TEC6_neg',
        'TEC2_neg+TEC4_neg+TEC5_neg',
        'TEC1_neg+TEC3_neg+TEC4_neg',
        'TEC6_neg+TEC7_neg+TEC9_neg',
        'TEC5_neg+TEC6_neg+TEC8_neg',
        'TEC4_neg+TEC5_neg+TEC7_neg',
        'TEC3_neg+TEC4_neg+TEC6_neg',
        'TEC2_neg+TEC3_neg+TEC5_neg',
        'TEC1_neg+TEC2_neg+TEC4_neg',
        'MTOB6+TEC1_pos+TEC2_pos',
        'MTOB6+TEC1_neg+TEC2_neg',
        'MTOB6+MTOB5+TEC1_pos',
        'MTOB6+MTOB5+TEC1_neg'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.simpleCosmicBONSeeds = cms.EDProducer("SimpleCosmicBONSeeder",
    ClusterChargeCheck = cms.PSet(
        Thresholds = cms.PSet(
            TEC = cms.int32(0),
            TIB = cms.int32(0),
            TID = cms.int32(0),
            TOB = cms.int32(0)
        ),
        checkCharge = cms.bool(False),
        matchedRecHitsUseAnd = cms.bool(True)
    ),
    ClusterCheckPSet = cms.PSet(
        ClusterCollectionLabel = cms.InputTag("siStripClusters"),
        DontCountDetsAboveNClusters = cms.uint32(20),
        MaxNumberOfCosmicClusters = cms.uint32(300),
        MaxNumberOfPixelClusters = cms.uint32(300),
        PixelClusterCollectionLabel = cms.InputTag("siPixelClusters"),
        doClusterCheck = cms.bool(True)
    ),
    HitsPerModuleCheck = cms.PSet(
        Thresholds = cms.PSet(
            TEC = cms.int32(20),
            TIB = cms.int32(20),
            TID = cms.int32(20),
            TOB = cms.int32(20)
        ),
        checkHitsPerModule = cms.bool(True)
    ),
    NegativeYOnly = cms.bool(False),
    PositiveYOnly = cms.bool(False),
    RegionPSet = cms.PSet(
        originHalfLength = cms.double(90.0),
        originRadius = cms.double(150.0),
        originZPosition = cms.double(0.0),
        pMin = cms.double(1.0),
        ptMin = cms.double(0.5)
    ),
    TTRHBuilder = cms.string('WithTrackAngle'),
    TripletsDebugLevel = cms.untracked.uint32(0),
    TripletsSrc = cms.InputTag("simpleCosmicBONSeedingLayers"),
    helixDebugLevel = cms.untracked.uint32(0),
    maxSeeds = cms.int32(20000),
    maxTriplets = cms.int32(50000),
    minimumGoodHitsInSeed = cms.int32(3),
    rescaleError = cms.double(1.0),
    seedDebugLevel = cms.untracked.uint32(0),
    seedOnMiddle = cms.bool(False),
    writeTriplets = cms.bool(False)
)


process.simpleCosmicBONSeedsBottom = cms.EDProducer("SimpleCosmicBONSeeder",
    ClusterChargeCheck = cms.PSet(
        Thresholds = cms.PSet(
            TEC = cms.int32(0),
            TIB = cms.int32(0),
            TID = cms.int32(0),
            TOB = cms.int32(0)
        ),
        checkCharge = cms.bool(False),
        matchedRecHitsUseAnd = cms.bool(True)
    ),
    ClusterCheckPSet = cms.PSet(
        ClusterCollectionLabel = cms.InputTag("siStripClustersBottom"),
        DontCountDetsAboveNClusters = cms.uint32(20),
        MaxNumberOfCosmicClusters = cms.uint32(150),
        MaxNumberOfPixelClusters = cms.uint32(300),
        PixelClusterCollectionLabel = cms.InputTag("siPixelClusters"),
        doClusterCheck = cms.bool(True)
    ),
    HitsPerModuleCheck = cms.PSet(
        Thresholds = cms.PSet(
            TEC = cms.int32(20),
            TIB = cms.int32(20),
            TID = cms.int32(20),
            TOB = cms.int32(20)
        ),
        checkHitsPerModule = cms.bool(True)
    ),
    NegativeYOnly = cms.bool(True),
    PositiveYOnly = cms.bool(False),
    RegionPSet = cms.PSet(
        originHalfLength = cms.double(90.0),
        originRadius = cms.double(150.0),
        originZPosition = cms.double(0.0),
        pMin = cms.double(1.0),
        ptMin = cms.double(0.5)
    ),
    TTRHBuilder = cms.string('WithTrackAngle'),
    TripletsDebugLevel = cms.untracked.uint32(0),
    TripletsSrc = cms.InputTag("simpleCosmicBONSeedingLayersBottom"),
    helixDebugLevel = cms.untracked.uint32(0),
    maxSeeds = cms.int32(20000),
    maxTriplets = cms.int32(50000),
    minimumGoodHitsInSeed = cms.int32(3),
    rescaleError = cms.double(1.0),
    seedDebugLevel = cms.untracked.uint32(0),
    seedOnMiddle = cms.bool(False),
    writeTriplets = cms.bool(False)
)


process.simpleCosmicBONSeedsTop = cms.EDProducer("SimpleCosmicBONSeeder",
    ClusterChargeCheck = cms.PSet(
        Thresholds = cms.PSet(
            TEC = cms.int32(0),
            TIB = cms.int32(0),
            TID = cms.int32(0),
            TOB = cms.int32(0)
        ),
        checkCharge = cms.bool(False),
        matchedRecHitsUseAnd = cms.bool(True)
    ),
    ClusterCheckPSet = cms.PSet(
        ClusterCollectionLabel = cms.InputTag("siStripClustersTop"),
        DontCountDetsAboveNClusters = cms.uint32(20),
        MaxNumberOfCosmicClusters = cms.uint32(150),
        MaxNumberOfPixelClusters = cms.uint32(300),
        PixelClusterCollectionLabel = cms.InputTag("siPixelClusters"),
        doClusterCheck = cms.bool(True)
    ),
    HitsPerModuleCheck = cms.PSet(
        Thresholds = cms.PSet(
            TEC = cms.int32(20),
            TIB = cms.int32(20),
            TID = cms.int32(20),
            TOB = cms.int32(20)
        ),
        checkHitsPerModule = cms.bool(True)
    ),
    NegativeYOnly = cms.bool(False),
    PositiveYOnly = cms.bool(True),
    RegionPSet = cms.PSet(
        originHalfLength = cms.double(90.0),
        originRadius = cms.double(150.0),
        originZPosition = cms.double(0.0),
        pMin = cms.double(1.0),
        ptMin = cms.double(0.5)
    ),
    TTRHBuilder = cms.string('WithTrackAngle'),
    TripletsDebugLevel = cms.untracked.uint32(0),
    TripletsSrc = cms.InputTag("simpleCosmicBONSeedingLayersTop"),
    helixDebugLevel = cms.untracked.uint32(0),
    maxSeeds = cms.int32(20000),
    maxTriplets = cms.int32(50000),
    minimumGoodHitsInSeed = cms.int32(3),
    rescaleError = cms.double(1.0),
    seedDebugLevel = cms.untracked.uint32(0),
    seedOnMiddle = cms.bool(False),
    writeTriplets = cms.bool(False)
)


process.splittedTracksP5 = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('cosmic'),
    Fitter = cms.string('RKFittingSmoother'),
    GeometricInnerState = cms.bool(False),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    Propagator = cms.string('RungeKuttaTrackerPropagator'),
    SimpleMagneticField = cms.string(''),
    TTRHBuilder = cms.string('WithTrackAngle'),
    TrajectoryInEvent = cms.bool(True),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("cosmicTrackSplitting"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(False)
)


process.stripPairElectronHitDoublets = cms.EDProducer("HitPairEDProducer",
    clusterCheck = cms.InputTag("trackerClusterCheck"),
    layerPairs = cms.vuint32(0),
    maxElement = cms.uint32(1000000),
    maxElementTotal = cms.uint32(12000000),
    mightGet = cms.optional.untracked.vstring,
    produceIntermediateHitDoublets = cms.bool(False),
    produceSeedingHitSets = cms.bool(True),
    seedingLayers = cms.InputTag("stripPairElectronSeedLayers"),
    trackingRegions = cms.InputTag("stripPairElectronTrackingRegions"),
    trackingRegionsSeedingLayers = cms.InputTag("")
)


process.stripPairElectronSeedLayers = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(

    ),
    FPix = cms.PSet(

    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        maxRing = cms.int32(2),
        minRing = cms.int32(1),
        skipClusters = cms.InputTag("tripletElectronClusterMask"),
        useRingSlector = cms.bool(True)
    ),
    TIB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        skipClusters = cms.InputTag("tripletElectronClusterMask")
    ),
    TID = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutNone')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        maxRing = cms.int32(2),
        minRing = cms.int32(1),
        skipClusters = cms.InputTag("tripletElectronClusterMask"),
        useRingSlector = cms.bool(True)
    ),
    TOB = cms.PSet(

    ),
    layerList = cms.vstring(
        'TIB1+TIB2',
        'TIB1+TID1_pos',
        'TIB1+TID1_neg',
        'TID2_pos+TID3_pos',
        'TID2_neg+TID3_neg',
        'TEC1_pos+TEC2_pos',
        'TEC2_pos+TEC3_pos',
        'TEC3_pos+TEC4_pos',
        'TEC3_pos+TEC5_pos',
        'TEC1_neg+TEC2_neg',
        'TEC2_neg+TEC3_neg',
        'TEC3_neg+TEC4_neg',
        'TEC3_neg+TEC5_neg'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.stripPairElectronSeeds = cms.EDProducer("SeedCreatorFromRegionConsecutiveHitsEDProducer",
    MinOneOverPtError = cms.double(1),
    OriginTransverseErrorMultiplier = cms.double(1),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('none')
    ),
    SeedMomentumForBOFF = cms.double(5),
    TTRHBuilder = cms.string('WithTrackAngle'),
    forceKinematicWithRegionDirection = cms.bool(False),
    magneticField = cms.string('ParabolicMf'),
    mightGet = cms.optional.untracked.vstring,
    propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    seedingHitSets = cms.InputTag("stripPairElectronHitDoublets")
)


process.stripPairElectronTrackingRegions = cms.EDProducer("GlobalTrackingRegionFromBeamSpotEDProducer",
    RegionPSet = cms.PSet(
        beamSpot = cms.InputTag("offlineBeamSpot"),
        nSigmaZ = cms.double(0),
        originHalfLength = cms.double(12.0),
        originRadius = cms.double(0.4),
        precise = cms.bool(True),
        ptMin = cms.double(1.0),
        useMultipleScattering = cms.bool(False)
    ),
    mightGet = cms.optional.untracked.vstring
)


process.tobTecStep = cms.EDProducer("TrackTfClassifier",
    beamspot = cms.InputTag("offlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mightGet = cms.optional.untracked.vstring,
    mva = cms.PSet(
        batchSize = cms.int32(16),
        tfDnnLabel = cms.string('trackSelectionTf')
    ),
    qualityCuts = cms.vdouble(-0.76, -0.65, -0.55),
    src = cms.InputTag("tobTecStepTracks"),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.tobTecStepClassifier1 = cms.EDProducer("TrackMVAClassifierDetached",
    beamspot = cms.InputTag("offlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mightGet = cms.optional.untracked.vstring,
    mva = cms.PSet(
        GBRForestFileName = cms.string(''),
        GBRForestLabel = cms.string('MVASelectorIter6_13TeV')
    ),
    qualityCuts = cms.vdouble(-0.6, -0.45, -0.3),
    src = cms.InputTag("tobTecStepTracks"),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.tobTecStepClassifier2 = cms.EDProducer("TrackMVAClassifierPrompt",
    beamspot = cms.InputTag("offlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mightGet = cms.optional.untracked.vstring,
    mva = cms.PSet(
        GBRForestFileName = cms.string(''),
        GBRForestLabel = cms.string('MVASelectorIter0_13TeV')
    ),
    qualityCuts = cms.vdouble(0.0, 0.0, 0.0),
    src = cms.InputTag("tobTecStepTracks"),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.tobTecStepClusters = cms.EDProducer("TrackClusterRemover",
    TrackQuality = cms.string('highPurity'),
    maxChi2 = cms.double(9.0),
    mightGet = cms.optional.untracked.vstring,
    minNumberOfLayersWithMeasBeforeFiltering = cms.int32(0),
    oldClusterRemovalInfo = cms.InputTag("pixelLessStepClusters"),
    overrideTrkQuals = cms.InputTag(""),
    pixelClusters = cms.InputTag("siPixelClusters"),
    stripClusters = cms.InputTag("siStripClusters"),
    trackClassifier = cms.InputTag("pixelLessStep","QualityMasks"),
    trajectories = cms.InputTag("pixelLessStepTracks")
)


process.tobTecStepHitDoubletsPair = cms.EDProducer("HitPairEDProducer",
    clusterCheck = cms.InputTag("trackerClusterCheck"),
    layerPairs = cms.vuint32(0),
    maxElement = cms.uint32(1000000),
    maxElementTotal = cms.uint32(12000000),
    mightGet = cms.optional.untracked.vstring,
    produceIntermediateHitDoublets = cms.bool(False),
    produceSeedingHitSets = cms.bool(True),
    seedingLayers = cms.InputTag("tobTecStepSeedLayersPair"),
    trackingRegions = cms.InputTag("tobTecStepTrackingRegionsPair"),
    trackingRegionsSeedingLayers = cms.InputTag("")
)


process.tobTecStepHitDoubletsTripl = cms.EDProducer("HitPairEDProducer",
    clusterCheck = cms.InputTag("trackerClusterCheck"),
    layerPairs = cms.vuint32(0),
    maxElement = cms.uint32(50000000),
    maxElementTotal = cms.uint32(50000000),
    mightGet = cms.optional.untracked.vstring,
    produceIntermediateHitDoublets = cms.bool(True),
    produceSeedingHitSets = cms.bool(False),
    seedingLayers = cms.InputTag("tobTecStepSeedLayersTripl"),
    trackingRegions = cms.InputTag("tobTecStepTrackingRegionsTripl"),
    trackingRegionsSeedingLayers = cms.InputTag("")
)


process.tobTecStepHitTripletsTripl = cms.EDProducer("MultiHitFromChi2EDProducer",
    ClusterShapeHitFilterName = cms.string('ClusterShapeHitFilter'),
    TTRHBuilder = cms.string('WithTrackAngle'),
    chi2VsPtCut = cms.bool(True),
    chi2_cuts = cms.vdouble(3, 4, 5, 5),
    detIdsToDebug = cms.vint32(0, 0, 0),
    doublets = cms.InputTag("tobTecStepHitDoubletsTripl"),
    extraHitRPhitolerance = cms.double(0),
    extraHitRZtolerance = cms.double(0),
    extraPhiKDBox = cms.double(0.01),
    extraRKDBox = cms.double(0.2),
    extraZKDBox = cms.double(0.2),
    fnSigmaRZ = cms.double(2),
    maxChi2 = cms.double(5),
    maxElement = cms.uint32(1000000),
    mightGet = cms.optional.untracked.vstring,
    phiPreFiltering = cms.double(0.3),
    pt_interv = cms.vdouble(0.4, 0.7, 1, 2),
    refitHits = cms.bool(True),
    useFixedPreFiltering = cms.bool(False)
)


process.tobTecStepMasks = cms.EDProducer("FastTrackerRecHitMaskProducer",
    TrackQuality = cms.string('highPurity'),
    minNumberOfLayersWithMeasBeforeFiltering = cms.int32(0),
    minNumberOfLayersWithMeasBeforeFiltering_ = cms.int32(0),
    oldHitRemovalInfo = cms.InputTag("pixelLessStepMasks"),
    recHits = cms.InputTag("fastTrackerRecHits"),
    trackClassifier = cms.InputTag("pixelLessStep","QualityMasks"),
    trackQuality = cms.string('hightPurity'),
    trajectories = cms.InputTag("pixelLessStepTracks")
)


process.tobTecStepSeedLayers = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(

    ),
    FPix = cms.PSet(

    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTiny')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        maxRing = cms.int32(5),
        minRing = cms.int32(5),
        skipClusters = cms.InputTag("tobTecStepClusters"),
        useRingSlector = cms.bool(True)
    ),
    TIB = cms.PSet(

    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTiny')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        skipClusters = cms.InputTag("tobTecStepClusters")
    ),
    layerList = cms.vstring(
        'TOB1+TOB2',
        'TOB1+TEC1_pos',
        'TOB1+TEC1_neg',
        'TEC1_pos+TEC2_pos',
        'TEC2_pos+TEC3_pos',
        'TEC3_pos+TEC4_pos',
        'TEC4_pos+TEC5_pos',
        'TEC5_pos+TEC6_pos',
        'TEC6_pos+TEC7_pos',
        'TEC1_neg+TEC2_neg',
        'TEC2_neg+TEC3_neg',
        'TEC3_neg+TEC4_neg',
        'TEC4_neg+TEC5_neg',
        'TEC5_neg+TEC6_neg',
        'TEC6_neg+TEC7_neg'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.tobTecStepSeedLayersPair = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(

    ),
    FPix = cms.PSet(

    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTight')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        maxRing = cms.int32(5),
        minRing = cms.int32(5),
        skipClusters = cms.InputTag("tobTecStepClusters"),
        useRingSlector = cms.bool(True)
    ),
    TIB = cms.PSet(

    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTight')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        skipClusters = cms.InputTag("tobTecStepClusters")
    ),
    layerList = cms.vstring(
        'TOB1+TEC1_pos',
        'TOB1+TEC1_neg',
        'TEC1_pos+TEC2_pos',
        'TEC1_neg+TEC2_neg',
        'TEC2_pos+TEC3_pos',
        'TEC2_neg+TEC3_neg',
        'TEC3_pos+TEC4_pos',
        'TEC3_neg+TEC4_neg',
        'TEC4_pos+TEC5_pos',
        'TEC4_neg+TEC5_neg',
        'TEC5_pos+TEC6_pos',
        'TEC5_neg+TEC6_neg',
        'TEC6_pos+TEC7_pos',
        'TEC6_neg+TEC7_neg'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.tobTecStepSeedLayersTripl = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(

    ),
    FPix = cms.PSet(

    ),
    MTEC = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTight')
        ),
        maxRing = cms.int32(7),
        minRing = cms.int32(6),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
        skipClusters = cms.InputTag("tobTecStepClusters"),
        useRingSlector = cms.bool(True)
    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTight')
        ),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
        skipClusters = cms.InputTag("tobTecStepClusters")
    ),
    TEC = cms.PSet(

    ),
    TIB = cms.PSet(

    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('SiStripClusterChargeCutTight')
        ),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        skipClusters = cms.InputTag("tobTecStepClusters")
    ),
    layerList = cms.vstring(
        'TOB1+TOB2+MTOB3',
        'TOB1+TOB2+MTOB4',
        'TOB1+TOB2+MTEC1_pos',
        'TOB1+TOB2+MTEC1_neg'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.tobTecStepSeeds = cms.EDProducer("SeedCombiner",
    seedCollections = cms.VInputTag("tobTecStepSeedsTripl", "tobTecStepSeedsPair")
)


process.tobTecStepSeedsPair = cms.EDProducer("SeedCreatorFromRegionConsecutiveHitsEDProducer",
    MinOneOverPtError = cms.double(1),
    OriginTransverseErrorMultiplier = cms.double(1),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('CombinedSeedComparitor'),
        comparitors = cms.VPSet(
            cms.PSet(
                ClusterShapeCacheSrc = cms.InputTag("siPixelClusterShapeCache"),
                ClusterShapeHitFilterName = cms.string('tobTecStepClusterShapeHitFilter'),
                ComponentName = cms.string('PixelClusterShapeSeedComparitor'),
                FilterAtHelixStage = cms.bool(True),
                FilterPixelHits = cms.bool(False),
                FilterStripHits = cms.bool(True)
            ),
            cms.PSet(
                ComponentName = cms.string('StripSubClusterShapeSeedFilter'),
                FilterAtHelixStage = cms.bool(False),
                label = cms.untracked.string('Seeds'),
                layerMask = cms.PSet(

                ),
                maxNSat = cms.uint32(3),
                maxTrimmedSizeDiffNeg = cms.double(1.0),
                maxTrimmedSizeDiffPos = cms.double(0.7),
                seedCutMIPs = cms.double(0.35),
                seedCutSN = cms.double(7.0),
                subclusterCutMIPs = cms.double(0.45),
                subclusterCutSN = cms.double(12.0),
                subclusterWindow = cms.double(0.7),
                trimMaxADC = cms.double(30.0),
                trimMaxFracNeigh = cms.double(0.25),
                trimMaxFracTotal = cms.double(0.15)
            )
        ),
        mode = cms.string('and')
    ),
    SeedMomentumForBOFF = cms.double(5),
    TTRHBuilder = cms.string('WithTrackAngle'),
    forceKinematicWithRegionDirection = cms.bool(False),
    magneticField = cms.string('ParabolicMf'),
    mightGet = cms.optional.untracked.vstring,
    propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    seedingHitSets = cms.InputTag("tobTecStepHitDoubletsPair")
)


process.tobTecStepSeedsTripl = cms.EDProducer("SeedCreatorFromRegionConsecutiveHitsEDProducer",
    MinOneOverPtError = cms.double(1),
    OriginTransverseErrorMultiplier = cms.double(1),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('CombinedSeedComparitor'),
        comparitors = cms.VPSet(
            cms.PSet(
                ClusterShapeCacheSrc = cms.InputTag("siPixelClusterShapeCache"),
                ClusterShapeHitFilterName = cms.string('tobTecStepClusterShapeHitFilter'),
                ComponentName = cms.string('PixelClusterShapeSeedComparitor'),
                FilterAtHelixStage = cms.bool(True),
                FilterPixelHits = cms.bool(False),
                FilterStripHits = cms.bool(True)
            ),
            cms.PSet(
                ComponentName = cms.string('StripSubClusterShapeSeedFilter'),
                FilterAtHelixStage = cms.bool(False),
                label = cms.untracked.string('Seeds'),
                layerMask = cms.PSet(

                ),
                maxNSat = cms.uint32(3),
                maxTrimmedSizeDiffNeg = cms.double(1.0),
                maxTrimmedSizeDiffPos = cms.double(0.7),
                seedCutMIPs = cms.double(0.35),
                seedCutSN = cms.double(7.0),
                subclusterCutMIPs = cms.double(0.45),
                subclusterCutSN = cms.double(12.0),
                subclusterWindow = cms.double(0.7),
                trimMaxADC = cms.double(30.0),
                trimMaxFracNeigh = cms.double(0.25),
                trimMaxFracTotal = cms.double(0.15)
            )
        ),
        mode = cms.string('and')
    ),
    SeedMomentumForBOFF = cms.double(5),
    TTRHBuilder = cms.string('WithTrackAngle'),
    forceKinematicWithRegionDirection = cms.bool(False),
    magneticField = cms.string('ParabolicMf'),
    mightGet = cms.optional.untracked.vstring,
    propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    seedingHitSets = cms.InputTag("tobTecStepHitTripletsTripl")
)


process.tobTecStepTrackCandidates = cms.EDProducer("CkfTrackCandidateMaker",
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    RedundantSeedCleaner = cms.string('CachingSeedCleanerBySharedInput'),
    TrajectoryBuilderPSet = cms.PSet(
        refToPSet_ = cms.string('tobTecStepTrajectoryBuilder')
    ),
    TrajectoryCleaner = cms.string('tobTecStepTrajectoryCleanerBySharedHits'),
    TransientInitialStateEstimatorParameters = cms.PSet(
        numberMeasurementsForFit = cms.int32(4),
        propagatorAlongTISE = cms.string('PropagatorWithMaterial'),
        propagatorOppositeTISE = cms.string('PropagatorWithMaterialOpposite')
    ),
    cleanTrajectoryAfterInOut = cms.bool(True),
    clustersToSkip = cms.InputTag("tobTecStepClusters"),
    doSeedingRegionRebuilding = cms.bool(True),
    maxNSeeds = cms.uint32(500000),
    maxSeedsBeforeCleaning = cms.uint32(5000),
    mightGet = cms.optional.untracked.vstring,
    numHitsForSeedCleaner = cms.int32(50),
    onlyPixelHitsForSeedCleaner = cms.bool(False),
    phase2clustersToSkip = cms.InputTag(""),
    reverseTrajectories = cms.bool(False),
    src = cms.InputTag("tobTecStepSeeds"),
    useHitsSplitting = cms.bool(True)
)


process.tobTecStepTrackCandidatesMkFit = cms.EDProducer("MkFitProducer",
    backwardFitInCMSSW = cms.bool(False),
    buildingRoutine = cms.string('cloneEngine'),
    clustersToSkip = cms.InputTag("tobTecStepClusters"),
    config = cms.ESInputTag("","tobTecStepTrackCandidatesMkFitConfig"),
    eventOfHits = cms.InputTag("mkFitEventOfHits"),
    limitConcurrency = cms.untracked.bool(False),
    mightGet = cms.optional.untracked.vstring,
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutLoose')
    ),
    mkFitSilent = cms.untracked.bool(True),
    pixelHits = cms.InputTag("mkFitSiPixelHits"),
    removeDuplicates = cms.bool(True),
    seedCleaning = cms.bool(True),
    seeds = cms.InputTag("tobTecStepTrackCandidatesMkFitSeeds"),
    stripHits = cms.InputTag("mkFitSiStripHits")
)


process.tobTecStepTrackCandidatesMkFitSeeds = cms.EDProducer("MkFitSeedConverter",
    maxNSeeds = cms.uint32(500000),
    mightGet = cms.optional.untracked.vstring,
    seeds = cms.InputTag("tobTecStepSeeds"),
    ttrhBuilder = cms.ESInputTag("","WithTrackAngle")
)


process.tobTecStepTrackingRegionsPair = cms.EDProducer("GlobalTrackingRegionFromBeamSpotEDProducer",
    RegionPSet = cms.PSet(
        beamSpot = cms.InputTag("offlineBeamSpot"),
        nSigmaZ = cms.double(0),
        originHalfLength = cms.double(30.0),
        originRadius = cms.double(6.0),
        precise = cms.bool(True),
        ptMin = cms.double(0.6),
        useMultipleScattering = cms.bool(False)
    ),
    mightGet = cms.optional.untracked.vstring
)


process.tobTecStepTrackingRegionsTripl = cms.EDProducer("GlobalTrackingRegionFromBeamSpotEDProducer",
    RegionPSet = cms.PSet(
        beamSpot = cms.InputTag("offlineBeamSpot"),
        nSigmaZ = cms.double(0),
        originHalfLength = cms.double(20.0),
        originRadius = cms.double(3.5),
        precise = cms.bool(True),
        ptMin = cms.double(0.55),
        useMultipleScattering = cms.bool(False)
    ),
    mightGet = cms.optional.untracked.vstring
)


process.tobTecStepTracks = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('tobTecStep'),
    Fitter = cms.string('tobTecFlexibleKFFittingSmoother'),
    GeometricInnerState = cms.bool(False),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    Propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    TTRHBuilder = cms.string('WithAngleAndTemplateWithoutProbQ'),
    TrajectoryInEvent = cms.bool(False),
    alias = cms.untracked.string('ctfWithMaterialTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("tobTecStepTrackCandidates"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(True)
)


process.topBottomClusterInfoProducer = cms.EDProducer("TopBottomClusterInfoProducer",
    pixelClustersNew = cms.InputTag("siPixelClustersTop"),
    pixelClustersOld = cms.InputTag("siPixelClusters"),
    pixelHitsNew = cms.InputTag("siPixelRecHitsTop"),
    pixelHitsOld = cms.InputTag("siPixelRecHits"),
    stripClustersNew = cms.InputTag("siStripClustersTop"),
    stripClustersOld = cms.InputTag("siStripClusters"),
    stripMonoHitsNew = cms.InputTag("siStripMatchedRecHitsTop","rphiRecHit"),
    stripMonoHitsOld = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
    stripStereoHitsNew = cms.InputTag("siStripMatchedRecHitsTop","stereoRecHit"),
    stripStereoHitsOld = cms.InputTag("siStripMatchedRecHits","stereoRecHit")
)


process.topBottomClusterInfoProducerBottom = cms.EDProducer("TopBottomClusterInfoProducer",
    pixelClustersNew = cms.InputTag("siPixelClustersBottom"),
    pixelClustersOld = cms.InputTag("siPixelClusters"),
    pixelHitsNew = cms.InputTag("siPixelRecHitsBottom"),
    pixelHitsOld = cms.InputTag("siPixelRecHits"),
    stripClustersNew = cms.InputTag("siStripClustersBottom"),
    stripClustersOld = cms.InputTag("siStripClusters"),
    stripMonoHitsNew = cms.InputTag("siStripMatchedRecHitsBottom","rphiRecHit"),
    stripMonoHitsOld = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
    stripStereoHitsNew = cms.InputTag("siStripMatchedRecHitsBottom","stereoRecHit"),
    stripStereoHitsOld = cms.InputTag("siStripMatchedRecHits","stereoRecHit")
)


process.topBottomClusterInfoProducerTop = cms.EDProducer("TopBottomClusterInfoProducer",
    pixelClustersNew = cms.InputTag("siPixelClustersTop"),
    pixelClustersOld = cms.InputTag("siPixelClusters"),
    pixelHitsNew = cms.InputTag("siPixelRecHitsTop"),
    pixelHitsOld = cms.InputTag("siPixelRecHits"),
    stripClustersNew = cms.InputTag("siStripClustersTop"),
    stripClustersOld = cms.InputTag("siStripClusters"),
    stripMonoHitsNew = cms.InputTag("siStripMatchedRecHitsTop","rphiRecHit"),
    stripMonoHitsOld = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
    stripStereoHitsNew = cms.InputTag("siStripMatchedRecHitsTop","stereoRecHit"),
    stripStereoHitsOld = cms.InputTag("siStripMatchedRecHits","stereoRecHit")
)


process.trackClusterRemover = cms.EDProducer("TrackClusterRemover",
    TrackQuality = cms.string('highPurity'),
    maxChi2 = cms.double(30),
    mightGet = cms.optional.untracked.vstring,
    minNumberOfLayersWithMeasBeforeFiltering = cms.int32(0),
    oldClusterRemovalInfo = cms.InputTag(""),
    overrideTrkQuals = cms.InputTag(""),
    pixelClusters = cms.InputTag("siPixelClusters"),
    stripClusters = cms.InputTag("siStripClusters"),
    trackClassifier = cms.InputTag("","QualityMasks"),
    trajectories = cms.InputTag("")
)


process.trackExtrapolator = cms.EDProducer("TrackExtrapolator",
    trackQuality = cms.string('goodIterative'),
    trackSrc = cms.InputTag("generalTracks")
)


process.trackRefsForJets = cms.EDProducer("ChargedRefCandidateProducer",
    particleType = cms.string('pi+'),
    src = cms.InputTag("trackWithVertexRefSelector")
)


process.trackTfClassifier = cms.EDProducer("TrackTfClassifier",
    beamspot = cms.InputTag("offlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mightGet = cms.optional.untracked.vstring,
    mva = cms.PSet(
        batchSize = cms.int32(16),
        tfDnnLabel = cms.string('trackSelectionTf')
    ),
    qualityCuts = cms.vdouble(-0.7, 0.1, 0.7),
    src = cms.InputTag(""),
    vertices = cms.InputTag("firstStepPrimaryVertices")
)


process.trackerClusterCheck = cms.EDProducer("ClusterCheckerEDProducer",
    ClusterCollectionLabel = cms.InputTag("siStripClusters"),
    MaxNumberOfCosmicClusters = cms.uint32(400000),
    MaxNumberOfPixelClusters = cms.uint32(40000),
    PixelClusterCollectionLabel = cms.InputTag("siPixelClusters"),
    cut = cms.string('strip < 400000 && pixel < 40000 && (strip < 50000 + 10*pixel) && (pixel < 5000 + 0.1*strip)'),
    doClusterCheck = cms.bool(False),
    mightGet = cms.optional.untracked.vstring,
    silentClusterCheck = cms.untracked.bool(False)
)


process.trackerClusterCheckPreSplitting = cms.EDProducer("ClusterCheckerEDProducer",
    ClusterCollectionLabel = cms.InputTag("siStripClusters"),
    MaxNumberOfCosmicClusters = cms.uint32(400000),
    MaxNumberOfPixelClusters = cms.uint32(40000),
    PixelClusterCollectionLabel = cms.InputTag("siPixelClustersPreSplitting"),
    cut = cms.string('strip < 400000 && pixel < 40000 && (strip < 50000 + 10*pixel) && (pixel < 5000 + 0.1*strip)'),
    doClusterCheck = cms.bool(False),
    mightGet = cms.optional.untracked.vstring,
    silentClusterCheck = cms.untracked.bool(False)
)


process.tripletElectronClusterMask = cms.EDProducer("SeedClusterRemover",
    Common = cms.PSet(
        maxChi2 = cms.double(9.0)
    ),
    oldClusterRemovalInfo = cms.InputTag("pixelLessStepSeedClusterMask"),
    pixelClusters = cms.InputTag("siPixelClusters"),
    stripClusters = cms.InputTag("siStripClusters"),
    trajectories = cms.InputTag("tripletElectronSeeds")
)


process.tripletElectronHitDoublets = cms.EDProducer("HitPairEDProducer",
    clusterCheck = cms.InputTag("trackerClusterCheck"),
    layerPairs = cms.vuint32(0),
    maxElement = cms.uint32(50000000),
    maxElementTotal = cms.uint32(50000000),
    mightGet = cms.optional.untracked.vstring,
    produceIntermediateHitDoublets = cms.bool(True),
    produceSeedingHitSets = cms.bool(False),
    seedingLayers = cms.InputTag("tripletElectronSeedLayers"),
    trackingRegions = cms.InputTag("tripletElectronTrackingRegions"),
    trackingRegionsSeedingLayers = cms.InputTag("")
)


process.tripletElectronHitTriplets = cms.EDProducer("PixelTripletHLTEDProducer",
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('none')
    ),
    doublets = cms.InputTag("tripletElectronHitDoublets"),
    extraHitRPhitolerance = cms.double(0.032),
    extraHitRZtolerance = cms.double(0.037),
    maxElement = cms.uint32(1000000),
    mightGet = cms.optional.untracked.vstring,
    phiPreFiltering = cms.double(0.3),
    produceIntermediateHitTriplets = cms.bool(False),
    produceSeedingHitSets = cms.bool(True),
    useBending = cms.bool(True),
    useFixedPreFiltering = cms.bool(False),
    useMultScattering = cms.bool(True)
)


process.tripletElectronSeedLayers = cms.EDProducer("SeedingLayersEDProducer",
    BPix = cms.PSet(
        HitProducer = cms.string('siPixelRecHits'),
        TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4PixelTriplets'),
        skipClusters = cms.InputTag("pixelLessStepSeedClusterMask")
    ),
    FPix = cms.PSet(
        HitProducer = cms.string('siPixelRecHits'),
        TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4PixelTriplets'),
        skipClusters = cms.InputTag("pixelLessStepSeedClusterMask")
    ),
    MTEC = cms.PSet(

    ),
    MTIB = cms.PSet(

    ),
    MTID = cms.PSet(

    ),
    MTOB = cms.PSet(

    ),
    TEC = cms.PSet(

    ),
    TIB = cms.PSet(

    ),
    TID = cms.PSet(

    ),
    TOB = cms.PSet(

    ),
    layerList = cms.vstring(
        'BPix1+BPix2+BPix3',
        'BPix2+BPix3+BPix4',
        'BPix1+BPix3+BPix4',
        'BPix1+BPix2+BPix4',
        'BPix2+BPix3+FPix1_pos',
        'BPix2+BPix3+FPix1_neg',
        'BPix1+BPix2+FPix1_pos',
        'BPix1+BPix2+FPix1_neg',
        'BPix1+BPix3+FPix1_pos',
        'BPix1+BPix3+FPix1_neg',
        'BPix2+FPix1_pos+FPix2_pos',
        'BPix2+FPix1_neg+FPix2_neg',
        'BPix1+FPix1_pos+FPix2_pos',
        'BPix1+FPix1_neg+FPix2_neg',
        'BPix1+BPix2+FPix2_pos',
        'BPix1+BPix2+FPix2_neg',
        'FPix1_pos+FPix2_pos+FPix3_pos',
        'FPix1_neg+FPix2_neg+FPix3_neg',
        'BPix1+FPix2_pos+FPix3_pos',
        'BPix1+FPix2_neg+FPix3_neg',
        'BPix1+FPix1_pos+FPix3_pos',
        'BPix1+FPix1_neg+FPix3_neg'
    ),
    mightGet = cms.optional.untracked.vstring
)


process.tripletElectronSeeds = cms.EDProducer("SeedCreatorFromRegionConsecutiveHitsEDProducer",
    MinOneOverPtError = cms.double(1),
    OriginTransverseErrorMultiplier = cms.double(1),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('none')
    ),
    SeedMomentumForBOFF = cms.double(5),
    TTRHBuilder = cms.string('WithTrackAngle'),
    forceKinematicWithRegionDirection = cms.bool(False),
    magneticField = cms.string('ParabolicMf'),
    mightGet = cms.optional.untracked.vstring,
    propagator = cms.string('PropagatorWithMaterialParabolicMf'),
    seedingHitSets = cms.InputTag("tripletElectronHitTriplets")
)


process.tripletElectronTrackingRegions = cms.EDProducer("GlobalTrackingRegionFromBeamSpotEDProducer",
    RegionPSet = cms.PSet(
        beamSpot = cms.InputTag("offlineBeamSpot"),
        nSigmaZ = cms.double(4),
        originHalfLength = cms.double(0),
        originRadius = cms.double(0.02),
        precise = cms.bool(True),
        ptMin = cms.double(1.0),
        useMultipleScattering = cms.bool(False)
    ),
    mightGet = cms.optional.untracked.vstring
)


process.siPixelClustersPreSplitting = SwitchProducerCUDA(
    cpu = cms.EDProducer("SiPixelClusterProducer",
        ChannelThreshold = cms.int32(10),
        ClusterMode = cms.string('PixelThresholdClusterizer'),
        ClusterThreshold = cms.int32(4000),
        ClusterThreshold_L1 = cms.int32(4000),
        DropDuplicates = cms.bool(True),
        ElectronPerADCGain = cms.double(135),
        MissCalibrate = cms.bool(True),
        Phase2Calibration = cms.bool(False),
        Phase2DigiBaseline = cms.double(1200),
        Phase2KinkADC = cms.int32(8),
        Phase2ReadoutMode = cms.int32(-1),
        SeedThreshold = cms.int32(1000),
        SplitClusters = cms.bool(False),
        VCaltoElectronGain = cms.int32(1),
        VCaltoElectronGain_L1 = cms.int32(1),
        VCaltoElectronOffset = cms.int32(0),
        VCaltoElectronOffset_L1 = cms.int32(0),
        maxNumberOfClusters = cms.int32(-1),
        mightGet = cms.optional.untracked.vstring,
        payloadType = cms.string('Offline'),
        src = cms.InputTag("siPixelDigis")
    )
)


process.siPixelRecHitsPreSplitting = SwitchProducerCUDA(
    cpu = cms.EDProducer("SiPixelRecHitConverter",
        CPE = cms.string('PixelCPEGeneric'),
        VerboseLevel = cms.untracked.int32(0),
        src = cms.InputTag("siPixelClustersPreSplitting")
    )
)


process.siPixelRecHitsPreSplittingSoA = SwitchProducerCUDA(
    cpu = cms.EDAlias(
        siPixelRecHitsPreSplittingCPU = cms.VPSet(
            cms.PSet(
                type = cms.string('pixelTopologyPhase1TrackingRecHitSoAHost')
            ),
            cms.PSet(
                type = cms.string('uintAsHostProduct')
            )
        )
    )
)


process.firstStepGoodPrimaryVertices = cms.EDFilter("PrimaryVertexObjectFilter",
    filterParams = cms.PSet(
        maxRho = cms.double(2.0),
        maxZ = cms.double(15.0),
        minNdof = cms.double(25.0)
    ),
    src = cms.InputTag("firstStepPrimaryVertices")
)


process.hltFilter = cms.EDFilter("HLTHighLevel",
    HLTPaths = cms.vstring('*'),
    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
    andOr = cms.bool(True),
    eventSetupPathsKey = cms.string(''),
    throw = cms.bool(False)
)


process.jetsForCoreTracking = cms.EDFilter("CandPtrSelector",
    cut = cms.string('pt > 100 && abs(eta) < 2.5'),
    filter = cms.bool(False),
    src = cms.InputTag("ak4CaloJetsForTrk")
)


process.jetsForCoreTrackingBarrel = cms.EDFilter("CandPtrSelector",
    cut = cms.string('pt > 100 && abs(eta) < 2.5'),
    filter = cms.bool(False),
    src = cms.InputTag("ak4CaloJetsForTrk")
)


process.jetsForCoreTrackingEndcap = cms.EDFilter("CandPtrSelector",
    cut = cms.string('pt > 100 && abs(eta) > 1.4 && abs(eta) < 2.5'),
    filter = cms.bool(False),
    src = cms.InputTag("ak4CaloJetsForTrk")
)


process.jetsForCoreTrackingPreSplitting = cms.EDFilter("CandPtrSelector",
    cut = cms.string('pt > 100 && abs(eta) < 2.5'),
    filter = cms.bool(False),
    src = cms.InputTag("ak4CaloJetsForTrkPreSplitting")
)


process.recoTauPileUpVertices = cms.EDFilter("RecoTauPileUpVertexSelector",
    filter = cms.bool(False),
    minTrackSumPt = cms.double(5),
    src = cms.InputTag("offlinePrimaryVertices")
)


process.fevt = cms.EDAnalyzer("RecHitAnalyzer",
    ak4PFJetCollection = cms.InputTag("ak4PFJets"),
    eleCollection = cms.InputTag("gedGsfElectrons"),
    gedPhotonCollection = cms.InputTag("gedPhotons"),
    genJetCollection = cms.InputTag("ak4GenJets"),
    genParticleCollection = cms.InputTag("genParticles"),
    granularityMultiEta = cms.int32(5),
    granularityMultiPhi = cms.int32(5),
    ipTagInfoCollection = cms.InputTag("pfImpactParameterTagInfos"),
    isBoostedTop = cms.bool(False),
    isDebug = cms.bool(False),
    isMC = cms.bool(True),
    isSignal = cms.bool(True),
    isW = cms.bool(False),
    jetTagCollection = cms.InputTag("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
    maxJetEta = cms.double(2.4),
    metCollection = cms.InputTag("pfMet"),
    minJetPt = cms.double(20.0),
    mode = cms.string('JetLevel'),
    nJets = cms.int32(-1),
    parameters = cms.PSet(
        dRMatch = cms.double(0.4),
        jetThreshold = cms.double(15),
        jeta = cms.vdouble(0.8, 1.3, 1.9, 2.5),
        jpar = cms.vdouble(1.39, 1.26, 1.21, 1.23, 1.28),
        pjpar = cms.vdouble(-0.2586, 0.6173)
    ),
    pfCollection = cms.InputTag("particleFlow"),
    recoJetsForBTagging = cms.InputTag("ak4PFJetsCHS"),
    reducedEBRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    reducedEERecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    reducedHBHERecHitCollection = cms.InputTag("hbhereco"),
    rhoLabel = cms.InputTag("fixedGridRhoAll"),
    secVertexCollection = cms.InputTag("inclusiveCandidateSecondaryVertices"),
    siPixelRecHitCollection = cms.InputTag("siPixelRecHits"),
    siStripMatchedRecHitCollection = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
    siStripRphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
    siStripStereoRecHits = cms.InputTag("siStripMatchedRecHits","stereoRecHit"),
    siStripUnmatchedRphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHitUnmatched"),
    siStripUnmatchedStereoRecHits = cms.InputTag("siStripMatchedRecHits","stereoRecHitUnmatched"),
    srcJetResPhi = cms.string('AK4PFchs_phi'),
    srcJetResPt = cms.string('AK4PFchs_pt'),
    srcJetSF = cms.string('AK4PFchs'),
    srcLeptons = cms.VInputTag("gedGsfElectrons", "muons", "gedPhotons"),
    srcPFCandidates = cms.InputTag("particleFlow"),
    srcPfJets = cms.InputTag("ak4PFJets"),
    task = cms.string('dijet_tau_massregression_unphysical'),
    tauCollection = cms.InputTag("hpsPFTauProducer"),
    tauDecayMode = cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding"),
    tauMVAIsolationRaw = cms.InputTag("hpsPFTauDiscriminationByIsolationMVArun2v1DBoldDMwLTraw"),
    trackCollection = cms.InputTag("generalTracks"),
    trackRecHitCollection = cms.InputTag("generalTracks"),
    transTrackBuilder = cms.ESInputTag("","TransientTrackBuilder"),
    triggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
    triggerSummaryTag = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
    vertexCollection = cms.InputTag("offlinePrimaryVertices"),
    z0PVCut = cms.double(0.1)
)


process.MessageLogger = cms.Service("MessageLogger",
    cerr = cms.untracked.PSet(
        FwkReport = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            reportEvery = cms.untracked.int32(1)
        ),
        FwkSummary = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            reportEvery = cms.untracked.int32(1)
        ),
        INFO = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        Root_NoDictionary = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000)
        ),
        enable = cms.untracked.bool(True),
        enableStatistics = cms.untracked.bool(False),
        lineLength = cms.optional.untracked.int32,
        noLineBreaks = cms.optional.untracked.bool,
        noTimeStamps = cms.untracked.bool(False),
        resetStatistics = cms.untracked.bool(False),
        statisticsThreshold = cms.untracked.string('WARNING'),
        threshold = cms.untracked.string('INFO'),
        allowAnyLabel_=cms.optional.untracked.PSetTemplate(
            limit = cms.optional.untracked.int32,
            reportEvery = cms.untracked.int32(1),
            timespan = cms.optional.untracked.int32
        )
    ),
    cout = cms.untracked.PSet(
        enable = cms.untracked.bool(False),
        enableStatistics = cms.untracked.bool(False),
        lineLength = cms.optional.untracked.int32,
        noLineBreaks = cms.optional.untracked.bool,
        noTimeStamps = cms.optional.untracked.bool,
        resetStatistics = cms.untracked.bool(False),
        statisticsThreshold = cms.optional.untracked.string,
        threshold = cms.optional.untracked.string,
        allowAnyLabel_=cms.optional.untracked.PSetTemplate(
            limit = cms.optional.untracked.int32,
            reportEvery = cms.untracked.int32(1),
            timespan = cms.optional.untracked.int32
        )
    ),
    debugModules = cms.untracked.vstring(),
    default = cms.untracked.PSet(
        limit = cms.optional.untracked.int32,
        lineLength = cms.untracked.int32(80),
        noLineBreaks = cms.untracked.bool(False),
        noTimeStamps = cms.untracked.bool(False),
        reportEvery = cms.untracked.int32(1),
        statisticsThreshold = cms.untracked.string('INFO'),
        threshold = cms.untracked.string('INFO'),
        timespan = cms.optional.untracked.int32,
        allowAnyLabel_=cms.optional.untracked.PSetTemplate(
            limit = cms.optional.untracked.int32,
            reportEvery = cms.untracked.int32(1),
            timespan = cms.optional.untracked.int32
        )
    ),
    files = cms.untracked.PSet(
        allowAnyLabel_=cms.optional.untracked.PSetTemplate(
            enableStatistics = cms.untracked.bool(False),
            extension = cms.optional.untracked.string,
            filename = cms.optional.untracked.string,
            lineLength = cms.optional.untracked.int32,
            noLineBreaks = cms.optional.untracked.bool,
            noTimeStamps = cms.optional.untracked.bool,
            output = cms.optional.untracked.string,
            resetStatistics = cms.untracked.bool(False),
            statisticsThreshold = cms.optional.untracked.string,
            threshold = cms.optional.untracked.string,
            allowAnyLabel_=cms.optional.untracked.PSetTemplate(
                limit = cms.optional.untracked.int32,
                reportEvery = cms.untracked.int32(1),
                timespan = cms.optional.untracked.int32
            )
        )
    ),
    suppressDebug = cms.untracked.vstring(),
    suppressFwkInfo = cms.untracked.vstring(),
    suppressInfo = cms.untracked.vstring(),
    suppressWarning = cms.untracked.vstring(),
    allowAnyLabel_=cms.optional.untracked.PSetTemplate(
        limit = cms.optional.untracked.int32,
        reportEvery = cms.untracked.int32(1),
        timespan = cms.optional.untracked.int32
    )
)


process.TFileService = cms.Service("TFileService",
    fileName = cms.string('output.root')
)


process.AnalyticalPropagator = cms.ESProducer("AnalyticalPropagatorESProducer",
    ComponentName = cms.string('AnalyticalPropagator'),
    MaxDPhi = cms.double(1.6),
    PropagationDirection = cms.string('alongMomentum')
)


process.AnalyticalPropagatorParabolicMF = cms.ESProducer("AnalyticalPropagatorESProducer",
    ComponentName = cms.string('AnalyticalPropagatorParabolicMf'),
    MaxDPhi = cms.double(1.6),
    PropagationDirection = cms.string('alongMomentum'),
    SimpleMagneticField = cms.string('ParabolicMf')
)


process.AnyDirectionAnalyticalPropagator = cms.ESProducer("AnalyticalPropagatorESProducer",
    ComponentName = cms.string('AnyDirectionAnalyticalPropagator'),
    MaxDPhi = cms.double(1.6),
    PropagationDirection = cms.string('anyDirection')
)


process.BeamHaloMPropagatorAlong = cms.ESProducer("PropagatorWithMaterialESProducer",
    ComponentName = cms.string('BeamHaloMPropagatorAlong'),
    Mass = cms.double(0.105),
    MaxDPhi = cms.double(10000),
    PropagationDirection = cms.string('alongMomentum'),
    SimpleMagneticField = cms.string(''),
    ptMin = cms.double(-1.0),
    useRungeKutta = cms.bool(True)
)


process.BeamHaloMPropagatorOpposite = cms.ESProducer("PropagatorWithMaterialESProducer",
    ComponentName = cms.string('BeamHaloMPropagatorOpposite'),
    Mass = cms.double(0.105),
    MaxDPhi = cms.double(10000),
    PropagationDirection = cms.string('oppositeToMomentum'),
    SimpleMagneticField = cms.string(''),
    ptMin = cms.double(-1.0),
    useRungeKutta = cms.bool(True)
)


process.BeamHaloPropagatorAlong = cms.ESProducer("BeamHaloPropagatorESProducer",
    ComponentName = cms.string('BeamHaloPropagatorAlong'),
    CrossingTrackerPropagator = cms.string('BeamHaloSHPropagatorAlong'),
    EndCapTrackerPropagator = cms.string('BeamHaloMPropagatorAlong'),
    PropagationDirection = cms.string('alongMomentum')
)


process.BeamHaloPropagatorAny = cms.ESProducer("BeamHaloPropagatorESProducer",
    ComponentName = cms.string('BeamHaloPropagatorAny'),
    CrossingTrackerPropagator = cms.string('BeamHaloSHPropagatorAny'),
    EndCapTrackerPropagator = cms.string('BeamHaloMPropagatorAlong'),
    PropagationDirection = cms.string('anyDirection')
)


process.BeamHaloPropagatorOpposite = cms.ESProducer("BeamHaloPropagatorESProducer",
    ComponentName = cms.string('BeamHaloPropagatorOpposite'),
    CrossingTrackerPropagator = cms.string('BeamHaloSHPropagatorOpposite'),
    EndCapTrackerPropagator = cms.string('BeamHaloMPropagatorOpposite'),
    PropagationDirection = cms.string('oppositeToMomentum')
)


process.BeamHaloSHPropagatorAlong = cms.ESProducer("SteppingHelixPropagatorESProducer",
    ApplyRadX0Correction = cms.bool(True),
    AssumeNoMaterial = cms.bool(False),
    ComponentName = cms.string('BeamHaloSHPropagatorAlong'),
    NoErrorPropagation = cms.bool(False),
    PropagationDirection = cms.string('alongMomentum'),
    SetVBFPointer = cms.bool(False),
    VBFName = cms.string('VolumeBasedMagneticField'),
    debug = cms.bool(False),
    endcapShiftInZNeg = cms.double(0.0),
    endcapShiftInZPos = cms.double(0.0),
    returnTangentPlane = cms.bool(True),
    sendLogWarning = cms.bool(False),
    useEndcapShiftsInZ = cms.bool(False),
    useInTeslaFromMagField = cms.bool(False),
    useIsYokeFlag = cms.bool(True),
    useMagVolumes = cms.bool(True),
    useMatVolumes = cms.bool(True),
    useTuningForL2Speed = cms.bool(False)
)


process.BeamHaloSHPropagatorAny = cms.ESProducer("SteppingHelixPropagatorESProducer",
    ApplyRadX0Correction = cms.bool(True),
    AssumeNoMaterial = cms.bool(False),
    ComponentName = cms.string('BeamHaloSHPropagatorAny'),
    NoErrorPropagation = cms.bool(False),
    PropagationDirection = cms.string('anyDirection'),
    SetVBFPointer = cms.bool(False),
    VBFName = cms.string('VolumeBasedMagneticField'),
    debug = cms.bool(False),
    endcapShiftInZNeg = cms.double(0.0),
    endcapShiftInZPos = cms.double(0.0),
    returnTangentPlane = cms.bool(True),
    sendLogWarning = cms.bool(False),
    useEndcapShiftsInZ = cms.bool(False),
    useInTeslaFromMagField = cms.bool(False),
    useIsYokeFlag = cms.bool(True),
    useMagVolumes = cms.bool(True),
    useMatVolumes = cms.bool(True),
    useTuningForL2Speed = cms.bool(False)
)


process.BeamHaloSHPropagatorOpposite = cms.ESProducer("SteppingHelixPropagatorESProducer",
    ApplyRadX0Correction = cms.bool(True),
    AssumeNoMaterial = cms.bool(False),
    ComponentName = cms.string('BeamHaloSHPropagatorOpposite'),
    NoErrorPropagation = cms.bool(False),
    PropagationDirection = cms.string('oppositeToMomentum'),
    SetVBFPointer = cms.bool(False),
    VBFName = cms.string('VolumeBasedMagneticField'),
    debug = cms.bool(False),
    endcapShiftInZNeg = cms.double(0.0),
    endcapShiftInZPos = cms.double(0.0),
    returnTangentPlane = cms.bool(True),
    sendLogWarning = cms.bool(False),
    useEndcapShiftsInZ = cms.bool(False),
    useInTeslaFromMagField = cms.bool(False),
    useIsYokeFlag = cms.bool(True),
    useMagVolumes = cms.bool(True),
    useMatVolumes = cms.bool(True),
    useTuningForL2Speed = cms.bool(False)
)


process.CSCGeometryESModule = cms.ESProducer("CSCGeometryESModule",
    alignmentsLabel = cms.string(''),
    appendToDataLabel = cms.string(''),
    applyAlignment = cms.bool(True),
    debugV = cms.untracked.bool(False),
    fromDD4hep = cms.bool(False),
    fromDDD = cms.bool(False),
    useCentreTIOffsets = cms.bool(False),
    useGangedStripsInME1a = cms.bool(False),
    useOnlyWiresInME1a = cms.bool(False),
    useRealWireGeometry = cms.bool(True)
)


process.CaloGeometryBuilder = cms.ESProducer("CaloGeometryBuilder",
    SelectedCalos = cms.vstring(
        'HCAL',
        'ZDC',
        'EcalBarrel',
        'EcalEndcap',
        'EcalPreshower',
        'TOWER'
    )
)


process.CaloTopologyBuilder = cms.ESProducer("CaloTopologyBuilder")


process.CaloTowerGeometryFromDBEP = cms.ESProducer("CaloTowerGeometryFromDBEP",
    applyAlignment = cms.bool(False)
)


process.CaloTowerTopologyEP = cms.ESProducer("CaloTowerTopologyEP")


process.CastorDbProducer = cms.ESProducer("CastorDbProducer",
    appendToDataLabel = cms.string('')
)


process.CastorGeometryFromDBEP = cms.ESProducer("CastorGeometryFromDBEP",
    applyAlignment = cms.bool(False)
)


process.Chi2MeasurementEstimator = cms.ESProducer("Chi2MeasurementEstimatorESProducer",
    ComponentName = cms.string('Chi2'),
    MaxChi2 = cms.double(30),
    MaxDisplacement = cms.double(0.5),
    MaxSagitta = cms.double(2),
    MinPtForHitRecoveryInGluedDet = cms.double(1000000000000),
    MinimalTolerance = cms.double(0.5),
    appendToDataLabel = cms.string(''),
    nSigma = cms.double(3)
)


process.Chi2MeasurementEstimatorForCDC = cms.ESProducer("Chi2MeasurementEstimatorESProducer",
    ComponentName = cms.string('Chi2MeasurementEstimatorForCDC'),
    MaxChi2 = cms.double(100.0),
    MaxDisplacement = cms.double(500),
    MaxSagitta = cms.double(-1),
    MinPtForHitRecoveryInGluedDet = cms.double(100000),
    MinimalTolerance = cms.double(0.5),
    appendToDataLabel = cms.string(''),
    nSigma = cms.double(4.0)
)


process.Chi2MeasurementEstimatorForP5 = cms.ESProducer("Chi2MeasurementEstimatorESProducer",
    ComponentName = cms.string('Chi2MeasurementEstimatorForP5'),
    MaxChi2 = cms.double(100.0),
    MaxDisplacement = cms.double(100),
    MaxSagitta = cms.double(-1),
    MinPtForHitRecoveryInGluedDet = cms.double(100000),
    MinimalTolerance = cms.double(0.5),
    appendToDataLabel = cms.string(''),
    nSigma = cms.double(4.0)
)


process.ClusterShapeHitFilterESProducer = cms.ESProducer("ClusterShapeHitFilterESProducer",
    ComponentName = cms.string('ClusterShapeHitFilter'),
    PixelShapeFile = cms.string('RecoPixelVertexing/PixelLowPtUtilities/data/pixelShapePhase1_noL1.par'),
    PixelShapeFileL1 = cms.string('RecoPixelVertexing/PixelLowPtUtilities/data/pixelShapePhase1_loose.par'),
    appendToDataLabel = cms.string(''),
    clusterChargeCut = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutNone')
    ),
    doPixelShapeCut = cms.bool(True),
    doStripShapeCut = cms.bool(True),
    isPhase2 = cms.bool(False)
)


process.DDCompactViewESProducer = cms.ESProducer("DDCompactViewESProducer",
    appendToDataLabel = cms.string('')
)


process.DDSpecParRegistryESProducer = cms.ESProducer("DDSpecParRegistryESProducer",
    appendToDataLabel = cms.string('')
)


process.DDVectorRegistryESProducer = cms.ESProducer("DDVectorRegistryESProducer",
    appendToDataLabel = cms.string('')
)


process.DTGeometryESModule = cms.ESProducer("DTGeometryESModule",
    DDDetector = cms.ESInputTag("",""),
    alignmentsLabel = cms.string(''),
    appendToDataLabel = cms.string(''),
    applyAlignment = cms.bool(True),
    attribute = cms.string('MuStructure'),
    fromDD4hep = cms.bool(False),
    fromDDD = cms.bool(False),
    value = cms.string('MuonBarrelDT')
)


process.DummyDetLayerGeometry = cms.ESProducer("DetLayerGeometryESProducer",
    ComponentName = cms.string('DummyDetLayerGeometry')
)


process.EcalBarrelGeometryFromDBEP = cms.ESProducer("EcalBarrelGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.EcalElectronicsMappingBuilder = cms.ESProducer("EcalElectronicsMappingBuilder")


process.EcalEndcapGeometryFromDBEP = cms.ESProducer("EcalEndcapGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.EcalLaserCorrectionService = cms.ESProducer("EcalLaserCorrectionService",
    maxExtrapolationTimeInSec = cms.uint32(0)
)


process.EcalLaserCorrectionServiceMC = cms.ESProducer("EcalLaserCorrectionServiceMC",
    appendToDataLabel = cms.string('')
)


process.EcalPreshowerGeometryFromDBEP = cms.ESProducer("EcalPreshowerGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.EcalTrigTowerConstituentsMapBuilder = cms.ESProducer("EcalTrigTowerConstituentsMapBuilder",
    MapFile = cms.untracked.string('Geometry/EcalMapping/data/EndCap_TTMap.txt')
)


process.FittingSmootherRKP5 = cms.ESProducer("KFFittingSmootherESProducer",
    BreakTrajWith2ConsecutiveMissing = cms.bool(False),
    ComponentName = cms.string('FittingSmootherRKP5'),
    EstimateCut = cms.double(20.0),
    Fitter = cms.string('RKFitter'),
    HighEtaSwitch = cms.double(5),
    LogPixelProbabilityCut = cms.double(0),
    MaxFractionOutliers = cms.double(0.3),
    MaxNumberOfOutliers = cms.int32(3),
    MinDof = cms.int32(2),
    MinNumberOfHits = cms.int32(4),
    MinNumberOfHitsHighEta = cms.int32(5),
    NoInvalidHitsBeginEnd = cms.bool(True),
    NoOutliersBeginEnd = cms.bool(False),
    RejectTracks = cms.bool(True),
    Smoother = cms.string('RKSmoother'),
    appendToDataLabel = cms.string('')
)


process.FlexibleKFFittingSmoother = cms.ESProducer("FlexibleKFFittingSmootherESProducer",
    ComponentName = cms.string('FlexibleKFFittingSmoother'),
    appendToDataLabel = cms.string(''),
    looperFitter = cms.string('LooperFittingSmoother'),
    standardFitter = cms.string('KFFittingSmootherWithOutliersRejectionAndRK')
)


process.GEMGeometryESModule = cms.ESProducer("GEMGeometryESModule",
    alignmentsLabel = cms.string(''),
    appendToDataLabel = cms.string(''),
    applyAlignment = cms.bool(True),
    fromDD4hep = cms.bool(False),
    fromDDD = cms.bool(False)
)


process.GlobalDetLayerGeometry = cms.ESProducer("GlobalDetLayerGeometryESProducer",
    ComponentName = cms.string('GlobalDetLayerGeometry')
)


process.GlobalTrackingGeometryESProducer = cms.ESProducer("GlobalTrackingGeometryESProducer")


process.HcalAlignmentEP = cms.ESProducer("HcalAlignmentEP")


process.HcalGeometryFromDBEP = cms.ESProducer("HcalGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.JetResolutionESProducer_AK4PFchs = cms.ESProducer("JetResolutionESProducer",
    label = cms.string('AK4PFchs')
)


process.JetResolutionESProducer_SF_AK4PFchs = cms.ESProducer("JetResolutionScaleFactorESProducer",
    label = cms.string('AK4PFchs')
)


process.KFFittingSmoother = cms.ESProducer("KFFittingSmootherESProducer",
    BreakTrajWith2ConsecutiveMissing = cms.bool(True),
    ComponentName = cms.string('KFFittingSmoother'),
    EstimateCut = cms.double(-1),
    Fitter = cms.string('KFFitter'),
    HighEtaSwitch = cms.double(5),
    LogPixelProbabilityCut = cms.double(0),
    MaxFractionOutliers = cms.double(0.3),
    MaxNumberOfOutliers = cms.int32(3),
    MinDof = cms.int32(2),
    MinNumberOfHits = cms.int32(5),
    MinNumberOfHitsHighEta = cms.int32(5),
    NoInvalidHitsBeginEnd = cms.bool(True),
    NoOutliersBeginEnd = cms.bool(False),
    RejectTracks = cms.bool(True),
    Smoother = cms.string('KFSmoother'),
    appendToDataLabel = cms.string('')
)


process.KFFittingSmootherBeamHalo = cms.ESProducer("KFFittingSmootherESProducer",
    BreakTrajWith2ConsecutiveMissing = cms.bool(True),
    ComponentName = cms.string('KFFittingSmootherBH'),
    EstimateCut = cms.double(-1),
    Fitter = cms.string('KFFitterBH'),
    HighEtaSwitch = cms.double(5),
    LogPixelProbabilityCut = cms.double(0),
    MaxFractionOutliers = cms.double(0.3),
    MaxNumberOfOutliers = cms.int32(3),
    MinDof = cms.int32(2),
    MinNumberOfHits = cms.int32(5),
    MinNumberOfHitsHighEta = cms.int32(5),
    NoInvalidHitsBeginEnd = cms.bool(True),
    NoOutliersBeginEnd = cms.bool(False),
    RejectTracks = cms.bool(True),
    Smoother = cms.string('KFSmootherBH'),
    appendToDataLabel = cms.string('')
)


process.KFFittingSmootherWithOutliersRejectionAndRK = cms.ESProducer("KFFittingSmootherESProducer",
    BreakTrajWith2ConsecutiveMissing = cms.bool(True),
    ComponentName = cms.string('KFFittingSmootherWithOutliersRejectionAndRK'),
    EstimateCut = cms.double(20.0),
    Fitter = cms.string('RKFitter'),
    HighEtaSwitch = cms.double(5),
    LogPixelProbabilityCut = cms.double(0),
    MaxFractionOutliers = cms.double(0.3),
    MaxNumberOfOutliers = cms.int32(3),
    MinDof = cms.int32(2),
    MinNumberOfHits = cms.int32(3),
    MinNumberOfHitsHighEta = cms.int32(5),
    NoInvalidHitsBeginEnd = cms.bool(True),
    NoOutliersBeginEnd = cms.bool(False),
    RejectTracks = cms.bool(True),
    Smoother = cms.string('RKSmoother'),
    appendToDataLabel = cms.string('')
)


process.KFSwitching1DUpdatorESProducer = cms.ESProducer("KFSwitching1DUpdatorESProducer",
    ComponentName = cms.string('KFSwitching1DUpdator'),
    doEndCap = cms.bool(False)
)


process.KFTrajectoryFitter = cms.ESProducer("KFTrajectoryFitterESProducer",
    ComponentName = cms.string('KFFitter'),
    Estimator = cms.string('Chi2'),
    Propagator = cms.string('PropagatorWithMaterial'),
    RecoGeometry = cms.string('GlobalDetLayerGeometry'),
    Updator = cms.string('KFUpdator'),
    appendToDataLabel = cms.string(''),
    minHits = cms.int32(3)
)


process.KFTrajectoryFitterBeamHalo = cms.ESProducer("KFTrajectoryFitterESProducer",
    ComponentName = cms.string('KFFitterBH'),
    Estimator = cms.string('Chi2'),
    Propagator = cms.string('BeamHaloPropagatorAlong'),
    RecoGeometry = cms.string('GlobalDetLayerGeometry'),
    Updator = cms.string('KFUpdator'),
    appendToDataLabel = cms.string(''),
    minHits = cms.int32(3)
)


process.KFTrajectorySmoother = cms.ESProducer("KFTrajectorySmootherESProducer",
    ComponentName = cms.string('KFSmoother'),
    Estimator = cms.string('Chi2'),
    Propagator = cms.string('PropagatorWithMaterial'),
    RecoGeometry = cms.string('GlobalDetLayerGeometry'),
    Updator = cms.string('KFUpdator'),
    appendToDataLabel = cms.string(''),
    errorRescaling = cms.double(100),
    minHits = cms.int32(3)
)


process.KFTrajectorySmootherBeamHalo = cms.ESProducer("KFTrajectorySmootherESProducer",
    ComponentName = cms.string('KFSmootherBH'),
    Estimator = cms.string('Chi2'),
    Propagator = cms.string('BeamHaloPropagatorAlong'),
    RecoGeometry = cms.string('GlobalDetLayerGeometry'),
    Updator = cms.string('KFUpdator'),
    appendToDataLabel = cms.string(''),
    errorRescaling = cms.double(100),
    minHits = cms.int32(3)
)


process.KFUpdatorESProducer = cms.ESProducer("KFUpdatorESProducer",
    ComponentName = cms.string('KFUpdator')
)


process.LooperFittingSmoother = cms.ESProducer("KFFittingSmootherESProducer",
    BreakTrajWith2ConsecutiveMissing = cms.bool(True),
    ComponentName = cms.string('LooperFittingSmoother'),
    EstimateCut = cms.double(20.0),
    Fitter = cms.string('LooperFitter'),
    HighEtaSwitch = cms.double(5),
    LogPixelProbabilityCut = cms.double(-14.0),
    MaxFractionOutliers = cms.double(0.3),
    MaxNumberOfOutliers = cms.int32(3),
    MinDof = cms.int32(2),
    MinNumberOfHits = cms.int32(3),
    MinNumberOfHitsHighEta = cms.int32(5),
    NoInvalidHitsBeginEnd = cms.bool(True),
    NoOutliersBeginEnd = cms.bool(False),
    RejectTracks = cms.bool(True),
    Smoother = cms.string('LooperSmoother'),
    appendToDataLabel = cms.string('')
)


process.LooperTrajectoryFitter = cms.ESProducer("KFTrajectoryFitterESProducer",
    ComponentName = cms.string('LooperFitter'),
    Estimator = cms.string('Chi2'),
    Propagator = cms.string('PropagatorWithMaterialForLoopers'),
    RecoGeometry = cms.string('GlobalDetLayerGeometry'),
    Updator = cms.string('KFUpdator'),
    appendToDataLabel = cms.string(''),
    minHits = cms.int32(3)
)


process.LooperTrajectorySmoother = cms.ESProducer("KFTrajectorySmootherESProducer",
    ComponentName = cms.string('LooperSmoother'),
    Estimator = cms.string('Chi2'),
    Propagator = cms.string('PropagatorWithMaterialForLoopers'),
    RecoGeometry = cms.string('GlobalDetLayerGeometry'),
    Updator = cms.string('KFUpdator'),
    appendToDataLabel = cms.string(''),
    errorRescaling = cms.double(10.0),
    minHits = cms.int32(3)
)


process.MRHChi2MeasurementEstimator = cms.ESProducer("MRHChi2MeasurementEstimatorESProducer",
    ComponentName = cms.string('MRHChi2'),
    MaxChi2 = cms.double(30.0),
    nSigma = cms.double(3.0)
)


process.MRHFittingSmoother = cms.ESProducer("KFFittingSmootherESProducer",
    BreakTrajWith2ConsecutiveMissing = cms.bool(True),
    ComponentName = cms.string('MRHFittingSmoother'),
    EstimateCut = cms.double(-1),
    Fitter = cms.string('MRHFitter'),
    HighEtaSwitch = cms.double(5),
    LogPixelProbabilityCut = cms.double(0),
    MaxFractionOutliers = cms.double(0.3),
    MaxNumberOfOutliers = cms.int32(3),
    MinDof = cms.int32(2),
    MinNumberOfHits = cms.int32(5),
    MinNumberOfHitsHighEta = cms.int32(5),
    NoInvalidHitsBeginEnd = cms.bool(True),
    NoOutliersBeginEnd = cms.bool(False),
    RejectTracks = cms.bool(True),
    Smoother = cms.string('MRHSmoother'),
    appendToDataLabel = cms.string('')
)


process.MRHTrajectoryFitter = cms.ESProducer("KFTrajectoryFitterESProducer",
    ComponentName = cms.string('MRHFitter'),
    Estimator = cms.string('MRHChi2'),
    Propagator = cms.string('RungeKuttaTrackerPropagator'),
    RecoGeometry = cms.string('GlobalDetLayerGeometry'),
    Updator = cms.string('KFUpdator'),
    appendToDataLabel = cms.string(''),
    minHits = cms.int32(3)
)


process.MRHTrajectorySmoother = cms.ESProducer("KFTrajectorySmootherESProducer",
    ComponentName = cms.string('MRHSmoother'),
    Estimator = cms.string('MRHChi2'),
    Propagator = cms.string('RungeKuttaTrackerPropagator'),
    RecoGeometry = cms.string('GlobalDetLayerGeometry'),
    Updator = cms.string('KFUpdator'),
    appendToDataLabel = cms.string(''),
    errorRescaling = cms.double(100),
    minHits = cms.int32(3)
)


process.MaterialPropagator = cms.ESProducer("PropagatorWithMaterialESProducer",
    ComponentName = cms.string('PropagatorWithMaterial'),
    Mass = cms.double(0.105),
    MaxDPhi = cms.double(1.6),
    PropagationDirection = cms.string('alongMomentum'),
    SimpleMagneticField = cms.string(''),
    ptMin = cms.double(-1.0),
    useRungeKutta = cms.bool(False)
)


process.MaterialPropagatorParabolicMF = cms.ESProducer("PropagatorWithMaterialESProducer",
    ComponentName = cms.string('PropagatorWithMaterialParabolicMf'),
    Mass = cms.double(0.105),
    MaxDPhi = cms.double(1.6),
    PropagationDirection = cms.string('alongMomentum'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    ptMin = cms.double(-1.0),
    useRungeKutta = cms.bool(False)
)


process.MeasurementTracker = cms.ESProducer("MeasurementTrackerESProducer",
    ComponentName = cms.string(''),
    DebugPixelModuleQualityDB = cms.untracked.bool(False),
    DebugPixelROCQualityDB = cms.untracked.bool(False),
    DebugStripAPVFiberQualityDB = cms.untracked.bool(False),
    DebugStripModuleQualityDB = cms.untracked.bool(False),
    DebugStripStripQualityDB = cms.untracked.bool(False),
    HitMatcher = cms.string('StandardMatcher'),
    MaskBadAPVFibers = cms.bool(True),
    Phase2StripCPE = cms.string(''),
    PixelCPE = cms.string('PixelCPEGeneric'),
    SiStripQualityLabel = cms.string(''),
    StripCPE = cms.string('StripCPEfromTrackAngle'),
    UsePixelModuleQualityDB = cms.bool(True),
    UsePixelROCQualityDB = cms.bool(True),
    UseStripAPVFiberQualityDB = cms.bool(True),
    UseStripModuleQualityDB = cms.bool(True),
    UseStripStripQualityDB = cms.bool(True),
    appendToDataLabel = cms.string(''),
    badStripCuts = cms.PSet(
        TEC = cms.PSet(
            maxBad = cms.uint32(4),
            maxConsecutiveBad = cms.uint32(2)
        ),
        TIB = cms.PSet(
            maxBad = cms.uint32(4),
            maxConsecutiveBad = cms.uint32(2)
        ),
        TID = cms.PSet(
            maxBad = cms.uint32(4),
            maxConsecutiveBad = cms.uint32(2)
        ),
        TOB = cms.PSet(
            maxBad = cms.uint32(4),
            maxConsecutiveBad = cms.uint32(2)
        )
    )
)


process.MeasurementTrackerBottom = cms.ESProducer("MeasurementTrackerESProducer",
    ComponentName = cms.string('MeasurementTrackerBottom'),
    DebugPixelModuleQualityDB = cms.untracked.bool(False),
    DebugPixelROCQualityDB = cms.untracked.bool(False),
    DebugStripAPVFiberQualityDB = cms.untracked.bool(False),
    DebugStripModuleQualityDB = cms.untracked.bool(False),
    DebugStripStripQualityDB = cms.untracked.bool(False),
    HitMatcher = cms.string('StandardMatcher'),
    MaskBadAPVFibers = cms.bool(True),
    Phase2StripCPE = cms.string(''),
    PixelCPE = cms.string('PixelCPEGeneric'),
    SiStripQualityLabel = cms.string(''),
    StripCPE = cms.string('StripCPEfromTrackAngle'),
    UsePixelModuleQualityDB = cms.bool(True),
    UsePixelROCQualityDB = cms.bool(True),
    UseStripAPVFiberQualityDB = cms.bool(True),
    UseStripModuleQualityDB = cms.bool(True),
    UseStripStripQualityDB = cms.bool(True),
    appendToDataLabel = cms.string(''),
    badStripCuts = cms.PSet(
        TEC = cms.PSet(
            maxBad = cms.uint32(4),
            maxConsecutiveBad = cms.uint32(2)
        ),
        TIB = cms.PSet(
            maxBad = cms.uint32(4),
            maxConsecutiveBad = cms.uint32(2)
        ),
        TID = cms.PSet(
            maxBad = cms.uint32(4),
            maxConsecutiveBad = cms.uint32(2)
        ),
        TOB = cms.PSet(
            maxBad = cms.uint32(4),
            maxConsecutiveBad = cms.uint32(2)
        )
    )
)


process.MeasurementTrackerTop = cms.ESProducer("MeasurementTrackerESProducer",
    ComponentName = cms.string('MeasurementTrackerTop'),
    DebugPixelModuleQualityDB = cms.untracked.bool(False),
    DebugPixelROCQualityDB = cms.untracked.bool(False),
    DebugStripAPVFiberQualityDB = cms.untracked.bool(False),
    DebugStripModuleQualityDB = cms.untracked.bool(False),
    DebugStripStripQualityDB = cms.untracked.bool(False),
    HitMatcher = cms.string('StandardMatcher'),
    MaskBadAPVFibers = cms.bool(True),
    Phase2StripCPE = cms.string(''),
    PixelCPE = cms.string('PixelCPEGeneric'),
    SiStripQualityLabel = cms.string(''),
    StripCPE = cms.string('StripCPEfromTrackAngle'),
    UsePixelModuleQualityDB = cms.bool(True),
    UsePixelROCQualityDB = cms.bool(True),
    UseStripAPVFiberQualityDB = cms.bool(True),
    UseStripModuleQualityDB = cms.bool(True),
    UseStripStripQualityDB = cms.bool(True),
    appendToDataLabel = cms.string(''),
    badStripCuts = cms.PSet(
        TEC = cms.PSet(
            maxBad = cms.uint32(4),
            maxConsecutiveBad = cms.uint32(2)
        ),
        TIB = cms.PSet(
            maxBad = cms.uint32(4),
            maxConsecutiveBad = cms.uint32(2)
        ),
        TID = cms.PSet(
            maxBad = cms.uint32(4),
            maxConsecutiveBad = cms.uint32(2)
        ),
        TOB = cms.PSet(
            maxBad = cms.uint32(4),
            maxConsecutiveBad = cms.uint32(2)
        )
    )
)


process.MuonDetLayerGeometryESProducer = cms.ESProducer("MuonDetLayerGeometryESProducer")


process.OppositeAnalyticalPropagator = cms.ESProducer("AnalyticalPropagatorESProducer",
    ComponentName = cms.string('AnalyticalPropagatorOpposite'),
    MaxDPhi = cms.double(1.6),
    PropagationDirection = cms.string('oppositeToMomentum')
)


process.OppositeAnalyticalPropagatorParabolicMF = cms.ESProducer("AnalyticalPropagatorESProducer",
    ComponentName = cms.string('AnalyticalPropagatorParabolicMfOpposite'),
    MaxDPhi = cms.double(1.6),
    PropagationDirection = cms.string('oppositeToMomentum'),
    SimpleMagneticField = cms.string('ParabolicMf')
)


process.OppositeMaterialPropagator = cms.ESProducer("PropagatorWithMaterialESProducer",
    ComponentName = cms.string('PropagatorWithMaterialOpposite'),
    Mass = cms.double(0.105),
    MaxDPhi = cms.double(1.6),
    PropagationDirection = cms.string('oppositeToMomentum'),
    SimpleMagneticField = cms.string(''),
    ptMin = cms.double(-1.0),
    useRungeKutta = cms.bool(False)
)


process.OppositeMaterialPropagatorParabolicMF = cms.ESProducer("PropagatorWithMaterialESProducer",
    ComponentName = cms.string('PropagatorWithMaterialParabolicMfOpposite'),
    Mass = cms.double(0.105),
    MaxDPhi = cms.double(1.6),
    PropagationDirection = cms.string('oppositeToMomentum'),
    SimpleMagneticField = cms.string('ParabolicMf'),
    ptMin = cms.double(-1.0),
    useRungeKutta = cms.bool(False)
)


process.ParabolicParametrizedMagneticFieldProducer = cms.ESProducer("AutoParametrizedMagneticFieldProducer",
    label = cms.untracked.string('ParabolicMf'),
    valueOverride = cms.int32(-1),
    version = cms.string('Parabolic')
)


process.PixelCPEGenericESProducer = cms.ESProducer("PixelCPEGenericESProducer",
    Alpha2Order = cms.bool(True),
    ClusterProbComputationFlag = cms.int32(0),
    ComponentName = cms.string('PixelCPEGeneric'),
    DoCosmics = cms.bool(False),
    EdgeClusterErrorX = cms.double(50),
    EdgeClusterErrorY = cms.double(85),
    IrradiationBiasCorrection = cms.bool(True),
    LoadTemplatesFromDB = cms.bool(True),
    MagneticFieldRecord = cms.ESInputTag("",""),
    NoTemplateErrorsWhenNoTrkAngles = cms.bool(False),
    SmallPitch = cms.bool(False),
    TruncatePixelCharge = cms.bool(True),
    UseErrorsFromTemplates = cms.bool(True),
    appendToDataLabel = cms.string(''),
    doLorentzFromAlignment = cms.bool(False),
    eff_charge_cut_highX = cms.double(1),
    eff_charge_cut_highY = cms.double(1),
    eff_charge_cut_lowX = cms.double(0),
    eff_charge_cut_lowY = cms.double(0),
    inflate_all_errors_no_trk_angle = cms.bool(False),
    inflate_errors = cms.bool(False),
    isPhase2 = cms.bool(False),
    lAOffset = cms.double(0),
    lAWidthBPix = cms.double(0),
    lAWidthFPix = cms.double(0),
    size_cutX = cms.double(3),
    size_cutY = cms.double(3),
    useLAFromDB = cms.bool(True),
    useLAWidthFromDB = cms.bool(True),
    xerr_barrel_l1 = cms.vdouble(0.00115, 0.0012, 0.00088),
    xerr_barrel_l1_def = cms.double(0.0103),
    xerr_barrel_ln = cms.vdouble(0.00115, 0.0012, 0.00088),
    xerr_barrel_ln_def = cms.double(0.0103),
    xerr_endcap = cms.vdouble(0.002, 0.002),
    xerr_endcap_def = cms.double(0.002),
    yerr_barrel_l1 = cms.vdouble(
        0.00375, 0.0023, 0.0025, 0.0025, 0.0023,
        0.0023, 0.0021, 0.0021, 0.0024
    ),
    yerr_barrel_l1_def = cms.double(0.0021),
    yerr_barrel_ln = cms.vdouble(
        0.00375, 0.0023, 0.0025, 0.0025, 0.0023,
        0.0023, 0.0021, 0.0021, 0.0024
    ),
    yerr_barrel_ln_def = cms.double(0.0021),
    yerr_endcap = cms.vdouble(0.0021),
    yerr_endcap_def = cms.double(0.00075)
)


process.PropagatorWithMaterialForLoopers = cms.ESProducer("PropagatorWithMaterialESProducer",
    ComponentName = cms.string('PropagatorWithMaterialForLoopers'),
    Mass = cms.double(0.1396),
    MaxDPhi = cms.double(4.0),
    PropagationDirection = cms.string('alongMomentum'),
    SimpleMagneticField = cms.string(''),
    ptMin = cms.double(-1),
    useOldAnalPropLogic = cms.bool(False),
    useRungeKutta = cms.bool(False)
)


process.PropagatorWithMaterialForLoopersOpposite = cms.ESProducer("PropagatorWithMaterialESProducer",
    ComponentName = cms.string('PropagatorWithMaterialForLoopersOpposite'),
    Mass = cms.double(0.1396),
    MaxDPhi = cms.double(4.0),
    PropagationDirection = cms.string('oppositeToMomentum'),
    SimpleMagneticField = cms.string(''),
    ptMin = cms.double(-1),
    useOldAnalPropLogic = cms.bool(False),
    useRungeKutta = cms.bool(False)
)


process.RK1DFittingSmoother = cms.ESProducer("KFFittingSmootherESProducer",
    BreakTrajWith2ConsecutiveMissing = cms.bool(True),
    ComponentName = cms.string('RK1DFittingSmoother'),
    EstimateCut = cms.double(-1),
    Fitter = cms.string('RK1DFitter'),
    HighEtaSwitch = cms.double(5),
    LogPixelProbabilityCut = cms.double(0),
    MaxFractionOutliers = cms.double(0.3),
    MaxNumberOfOutliers = cms.int32(3),
    MinDof = cms.int32(2),
    MinNumberOfHits = cms.int32(5),
    MinNumberOfHitsHighEta = cms.int32(5),
    NoInvalidHitsBeginEnd = cms.bool(True),
    NoOutliersBeginEnd = cms.bool(False),
    RejectTracks = cms.bool(True),
    Smoother = cms.string('RK1DSmoother'),
    appendToDataLabel = cms.string('')
)


process.RK1DTrajectoryFitter = cms.ESProducer("KFTrajectoryFitterESProducer",
    ComponentName = cms.string('RK1DFitter'),
    Estimator = cms.string('Chi2'),
    Propagator = cms.string('RungeKuttaTrackerPropagator'),
    RecoGeometry = cms.string('GlobalDetLayerGeometry'),
    Updator = cms.string('KFSwitching1DUpdator'),
    appendToDataLabel = cms.string(''),
    minHits = cms.int32(3)
)


process.RK1DTrajectorySmoother = cms.ESProducer("KFTrajectorySmootherESProducer",
    ComponentName = cms.string('RK1DSmoother'),
    Estimator = cms.string('Chi2'),
    Propagator = cms.string('RungeKuttaTrackerPropagator'),
    RecoGeometry = cms.string('GlobalDetLayerGeometry'),
    Updator = cms.string('KFSwitching1DUpdator'),
    appendToDataLabel = cms.string(''),
    errorRescaling = cms.double(100),
    minHits = cms.int32(3)
)


process.RKFittingSmoother = cms.ESProducer("KFFittingSmootherESProducer",
    BreakTrajWith2ConsecutiveMissing = cms.bool(True),
    ComponentName = cms.string('RKFittingSmoother'),
    EstimateCut = cms.double(-1),
    Fitter = cms.string('RKFitter'),
    HighEtaSwitch = cms.double(5),
    LogPixelProbabilityCut = cms.double(0),
    MaxFractionOutliers = cms.double(0.3),
    MaxNumberOfOutliers = cms.int32(3),
    MinDof = cms.int32(2),
    MinNumberOfHits = cms.int32(3),
    MinNumberOfHitsHighEta = cms.int32(5),
    NoInvalidHitsBeginEnd = cms.bool(True),
    NoOutliersBeginEnd = cms.bool(False),
    RejectTracks = cms.bool(True),
    Smoother = cms.string('RKSmoother'),
    appendToDataLabel = cms.string('')
)


process.RKOutliers1DFittingSmoother = cms.ESProducer("KFFittingSmootherESProducer",
    BreakTrajWith2ConsecutiveMissing = cms.bool(True),
    ComponentName = cms.string('RKOutliers1DFittingSmoother'),
    EstimateCut = cms.double(20.0),
    Fitter = cms.string('RK1DFitter'),
    HighEtaSwitch = cms.double(5),
    LogPixelProbabilityCut = cms.double(0),
    MaxFractionOutliers = cms.double(0.3),
    MaxNumberOfOutliers = cms.int32(3),
    MinDof = cms.int32(2),
    MinNumberOfHits = cms.int32(3),
    MinNumberOfHitsHighEta = cms.int32(5),
    NoInvalidHitsBeginEnd = cms.bool(True),
    NoOutliersBeginEnd = cms.bool(False),
    RejectTracks = cms.bool(True),
    Smoother = cms.string('RK1DSmoother'),
    appendToDataLabel = cms.string('')
)


process.RKTrajectoryFitter = cms.ESProducer("KFTrajectoryFitterESProducer",
    ComponentName = cms.string('RKFitter'),
    Estimator = cms.string('Chi2'),
    Propagator = cms.string('RungeKuttaTrackerPropagator'),
    RecoGeometry = cms.string('GlobalDetLayerGeometry'),
    Updator = cms.string('KFUpdator'),
    appendToDataLabel = cms.string(''),
    minHits = cms.int32(3)
)


process.RKTrajectorySmoother = cms.ESProducer("KFTrajectorySmootherESProducer",
    ComponentName = cms.string('RKSmoother'),
    Estimator = cms.string('Chi2'),
    Propagator = cms.string('RungeKuttaTrackerPropagator'),
    RecoGeometry = cms.string('GlobalDetLayerGeometry'),
    Updator = cms.string('KFUpdator'),
    appendToDataLabel = cms.string(''),
    errorRescaling = cms.double(100),
    minHits = cms.int32(3)
)


process.RPCGeometryESModule = cms.ESProducer("RPCGeometryESModule",
    fromDD4hep = cms.untracked.bool(False),
    fromDDD = cms.untracked.bool(False)
)


process.RungeKuttaTrackerPropagator = cms.ESProducer("PropagatorWithMaterialESProducer",
    ComponentName = cms.string('RungeKuttaTrackerPropagator'),
    Mass = cms.double(0.105),
    MaxDPhi = cms.double(1.6),
    PropagationDirection = cms.string('alongMomentum'),
    SimpleMagneticField = cms.string(''),
    ptMin = cms.double(-1.0),
    useRungeKutta = cms.bool(True)
)


process.SiPixelTemplateStoreESProducer = cms.ESProducer("SiPixelTemplateStoreESProducer",
    appendToDataLabel = cms.string('')
)


process.SiStripClusterizerConditionsESProducer = cms.ESProducer("SiStripClusterizerConditionsESProducer",
    Label = cms.string(''),
    QualityLabel = cms.string(''),
    appendToDataLabel = cms.string('')
)


process.SiStripRecHitMatcherESProducer = cms.ESProducer("SiStripRecHitMatcherESProducer",
    ComponentName = cms.string('StandardMatcher'),
    NSigmaInside = cms.double(3.0),
    PreFilter = cms.bool(False)
)


process.SteppingHelixPropagatorAlong = cms.ESProducer("SteppingHelixPropagatorESProducer",
    ApplyRadX0Correction = cms.bool(True),
    AssumeNoMaterial = cms.bool(False),
    ComponentName = cms.string('SteppingHelixPropagatorAlong'),
    NoErrorPropagation = cms.bool(False),
    PropagationDirection = cms.string('alongMomentum'),
    SetVBFPointer = cms.bool(False),
    VBFName = cms.string('VolumeBasedMagneticField'),
    debug = cms.bool(False),
    endcapShiftInZNeg = cms.double(0.0),
    endcapShiftInZPos = cms.double(0.0),
    returnTangentPlane = cms.bool(True),
    sendLogWarning = cms.bool(False),
    useEndcapShiftsInZ = cms.bool(False),
    useInTeslaFromMagField = cms.bool(False),
    useIsYokeFlag = cms.bool(True),
    useMagVolumes = cms.bool(True),
    useMatVolumes = cms.bool(True),
    useTuningForL2Speed = cms.bool(False)
)


process.SteppingHelixPropagatorAny = cms.ESProducer("SteppingHelixPropagatorESProducer",
    ApplyRadX0Correction = cms.bool(True),
    AssumeNoMaterial = cms.bool(False),
    ComponentName = cms.string('SteppingHelixPropagatorAny'),
    NoErrorPropagation = cms.bool(False),
    PropagationDirection = cms.string('anyDirection'),
    SetVBFPointer = cms.bool(False),
    VBFName = cms.string('VolumeBasedMagneticField'),
    debug = cms.bool(False),
    endcapShiftInZNeg = cms.double(0.0),
    endcapShiftInZPos = cms.double(0.0),
    returnTangentPlane = cms.bool(True),
    sendLogWarning = cms.bool(False),
    useEndcapShiftsInZ = cms.bool(False),
    useInTeslaFromMagField = cms.bool(False),
    useIsYokeFlag = cms.bool(True),
    useMagVolumes = cms.bool(True),
    useMatVolumes = cms.bool(True),
    useTuningForL2Speed = cms.bool(False)
)


process.SteppingHelixPropagatorOpposite = cms.ESProducer("SteppingHelixPropagatorESProducer",
    ApplyRadX0Correction = cms.bool(True),
    AssumeNoMaterial = cms.bool(False),
    ComponentName = cms.string('SteppingHelixPropagatorOpposite'),
    NoErrorPropagation = cms.bool(False),
    PropagationDirection = cms.string('oppositeToMomentum'),
    SetVBFPointer = cms.bool(False),
    VBFName = cms.string('VolumeBasedMagneticField'),
    debug = cms.bool(False),
    endcapShiftInZNeg = cms.double(0.0),
    endcapShiftInZPos = cms.double(0.0),
    returnTangentPlane = cms.bool(True),
    sendLogWarning = cms.bool(False),
    useEndcapShiftsInZ = cms.bool(False),
    useInTeslaFromMagField = cms.bool(False),
    useIsYokeFlag = cms.bool(True),
    useMagVolumes = cms.bool(True),
    useMatVolumes = cms.bool(True),
    useTuningForL2Speed = cms.bool(False)
)


process.StripCPEESProducer = cms.ESProducer("StripCPEESProducer",
    ComponentName = cms.string('SimpleStripCPE'),
    ComponentType = cms.string('SimpleStripCPE'),
    parameters = cms.PSet(

    )
)


process.StripCPEfromTrackAngleESProducer = cms.ESProducer("StripCPEESProducer",
    ComponentName = cms.string('StripCPEfromTrackAngle'),
    ComponentType = cms.string('StripCPEfromTrackAngle'),
    parameters = cms.PSet(
        mLC_P0 = cms.double(-0.326),
        mLC_P1 = cms.double(0.618),
        mLC_P2 = cms.double(0.3),
        mTEC_P0 = cms.double(-1.885),
        mTEC_P1 = cms.double(0.471),
        mTIB_P0 = cms.double(-0.742),
        mTIB_P1 = cms.double(0.202),
        mTID_P0 = cms.double(-1.427),
        mTID_P1 = cms.double(0.433),
        mTOB_P0 = cms.double(-1.026),
        mTOB_P1 = cms.double(0.253),
        maxChgOneMIP = cms.double(6000.0),
        useLegacyError = cms.bool(False)
    )
)


process.TTRHBuilderAngleAndTemplate = cms.ESProducer("TkTransientTrackingRecHitBuilderESProducer",
    ComponentName = cms.string('WithAngleAndTemplate'),
    ComputeCoarseLocalPositionFromDisk = cms.bool(False),
    Matcher = cms.string('StandardMatcher'),
    Phase2StripCPE = cms.string(''),
    PixelCPE = cms.string('PixelCPEClusterRepair'),
    StripCPE = cms.string('StripCPEfromTrackAngle'),
    appendToDataLabel = cms.string('')
)


process.TTRHBuilderAngleAndTemplateWithoutProbQ = cms.ESProducer("TkTransientTrackingRecHitBuilderESProducer",
    ComponentName = cms.string('WithAngleAndTemplateWithoutProbQ'),
    ComputeCoarseLocalPositionFromDisk = cms.bool(False),
    Matcher = cms.string('StandardMatcher'),
    Phase2StripCPE = cms.string(''),
    PixelCPE = cms.string('PixelCPEClusterRepairWithoutProbQ'),
    StripCPE = cms.string('StripCPEfromTrackAngle'),
    appendToDataLabel = cms.string('')
)


process.TTRHBuilderFast = cms.ESProducer("TkTransientTrackingRecHitBuilderESProducer",
    ComponentName = cms.string('WithoutAngleFast'),
    ComputeCoarseLocalPositionFromDisk = cms.bool(False),
    Matcher = cms.string('StandardMatcher'),
    Phase2StripCPE = cms.string(''),
    PixelCPE = cms.string('PixelCPEFast'),
    StripCPE = cms.string('StripCPEfromTrackAngle'),
    appendToDataLabel = cms.string('')
)


process.TrackerRecoGeometryESProducer = cms.ESProducer("TrackerRecoGeometryESProducer",
    usePhase2Stacks = cms.bool(False)
)


process.TransientTrackBuilderESProducer = cms.ESProducer("TransientTrackBuilderESProducer",
    ComponentName = cms.string('TransientTrackBuilder')
)


process.VolumeBasedMagneticFieldESProducer = cms.ESProducer("DD4hep_VolumeBasedMagneticFieldESProducerFromDB",
    debugBuilder = cms.untracked.bool(False),
    label = cms.untracked.string(''),
    valueOverride = cms.int32(-1)
)


process.ZdcGeometryFromDBEP = cms.ESProducer("ZdcGeometryFromDBEP",
    applyAlignment = cms.bool(False)
)


process.beamHaloNavigationSchoolESProducer = cms.ESProducer("NavigationSchoolESProducer",
    ComponentName = cms.string('BeamHaloNavigationSchool'),
    SimpleMagneticField = cms.string('')
)


process.caloSimulationParameters = cms.ESProducer("CaloSimParametersESModule",
    appendToDataLabel = cms.string(''),
    fromDD4hep = cms.bool(True)
)


process.clusterShapeHitFilterESProducer = cms.ESProducer("ClusterShapeHitFilterESProducer",
    ComponentName = cms.string(''),
    PixelShapeFile = cms.string('RecoPixelVertexing/PixelLowPtUtilities/data/pixelShapePhase0.par'),
    PixelShapeFileL1 = cms.string('RecoPixelVertexing/PixelLowPtUtilities/data/pixelShapePhase0.par'),
    appendToDataLabel = cms.string(''),
    clusterChargeCut = cms.PSet(
        value = cms.double(-1)
    ),
    doPixelShapeCut = cms.bool(True),
    doStripShapeCut = cms.bool(True),
    isPhase2 = cms.bool(False)
)


process.conv2StepFitterSmoother = cms.ESProducer("KFFittingSmootherESProducer",
    BreakTrajWith2ConsecutiveMissing = cms.bool(True),
    ComponentName = cms.string('conv2StepFitterSmoother'),
    EstimateCut = cms.double(30),
    Fitter = cms.string('RKFitter'),
    HighEtaSwitch = cms.double(5),
    LogPixelProbabilityCut = cms.double(0),
    MaxFractionOutliers = cms.double(0.3),
    MaxNumberOfOutliers = cms.int32(3),
    MinDof = cms.int32(2),
    MinNumberOfHits = cms.int32(3),
    MinNumberOfHitsHighEta = cms.int32(5),
    NoInvalidHitsBeginEnd = cms.bool(True),
    NoOutliersBeginEnd = cms.bool(False),
    RejectTracks = cms.bool(True),
    Smoother = cms.string('conv2StepRKSmoother'),
    appendToDataLabel = cms.string('')
)


process.conv2StepRKTrajectorySmoother = cms.ESProducer("KFTrajectorySmootherESProducer",
    ComponentName = cms.string('conv2StepRKSmoother'),
    Estimator = cms.string('Chi2'),
    Propagator = cms.string('RungeKuttaTrackerPropagator'),
    RecoGeometry = cms.string('GlobalDetLayerGeometry'),
    Updator = cms.string('KFUpdator'),
    appendToDataLabel = cms.string(''),
    errorRescaling = cms.double(10.0),
    minHits = cms.int32(3)
)


process.convStepChi2Est = cms.ESProducer("Chi2ChargeMeasurementEstimatorESProducer",
    ComponentName = cms.string('convStepChi2Est'),
    MaxChi2 = cms.double(30),
    MaxDisplacement = cms.double(100),
    MaxSagitta = cms.double(-1.0),
    MinPtForHitRecoveryInGluedDet = cms.double(1000000000000),
    MinimalTolerance = cms.double(0.5),
    appendToDataLabel = cms.string(''),
    clusterChargeCut = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutTight')
    ),
    nSigma = cms.double(3),
    pTChargeCutThreshold = cms.double(-1)
)


process.convStepFitterSmoother = cms.ESProducer("KFFittingSmootherESProducer",
    BreakTrajWith2ConsecutiveMissing = cms.bool(True),
    ComponentName = cms.string('convStepFitterSmoother'),
    EstimateCut = cms.double(30),
    Fitter = cms.string('RKFitter'),
    HighEtaSwitch = cms.double(5),
    LogPixelProbabilityCut = cms.double(0),
    MaxFractionOutliers = cms.double(0.3),
    MaxNumberOfOutliers = cms.int32(3),
    MinDof = cms.int32(2),
    MinNumberOfHits = cms.int32(3),
    MinNumberOfHitsHighEta = cms.int32(5),
    NoInvalidHitsBeginEnd = cms.bool(True),
    NoOutliersBeginEnd = cms.bool(False),
    RejectTracks = cms.bool(True),
    Smoother = cms.string('convStepRKSmoother'),
    appendToDataLabel = cms.string('')
)


process.convStepRKTrajectorySmoother = cms.ESProducer("KFTrajectorySmootherESProducer",
    ComponentName = cms.string('convStepRKSmoother'),
    Estimator = cms.string('Chi2'),
    Propagator = cms.string('RungeKuttaTrackerPropagator'),
    RecoGeometry = cms.string('GlobalDetLayerGeometry'),
    Updator = cms.string('KFUpdator'),
    appendToDataLabel = cms.string(''),
    errorRescaling = cms.double(10.0),
    minHits = cms.int32(3)
)


process.cosmicsNavigationSchoolESProducer = cms.ESProducer("SkippingLayerCosmicNavigationSchoolESProducer",
    ComponentName = cms.string('CosmicNavigationSchool'),
    allSelf = cms.bool(True),
    noPXB = cms.bool(False),
    noPXF = cms.bool(False),
    noTEC = cms.bool(False),
    noTIB = cms.bool(False),
    noTID = cms.bool(False),
    noTOB = cms.bool(False),
    selfSearch = cms.bool(True)
)


process.ctppsBeamParametersFromLHCInfoESSource = cms.ESProducer("CTPPSBeamParametersFromLHCInfoESSource",
    appendToDataLabel = cms.string(''),
    beamDivX45 = cms.double(0.1),
    beamDivX56 = cms.double(0.1),
    beamDivY45 = cms.double(0.1),
    beamDivY56 = cms.double(0.1),
    lhcInfoLabel = cms.string(''),
    vtxOffsetX45 = cms.double(0.01),
    vtxOffsetX56 = cms.double(0.01),
    vtxOffsetY45 = cms.double(0.01),
    vtxOffsetY56 = cms.double(0.01),
    vtxOffsetZ45 = cms.double(0.01),
    vtxOffsetZ56 = cms.double(0.01),
    vtxStddevX = cms.double(0.02),
    vtxStddevY = cms.double(0.02),
    vtxStddevZ = cms.double(0.02)
)


process.ctppsInterpolatedOpticalFunctionsESSource = cms.ESProducer("CTPPSInterpolatedOpticalFunctionsESSource",
    appendToDataLabel = cms.string(''),
    lhcInfoLabel = cms.string(''),
    opticsLabel = cms.string('')
)


process.detachedQuadStepChi2Est = cms.ESProducer("Chi2ChargeMeasurementEstimatorESProducer",
    ComponentName = cms.string('detachedQuadStepChi2Est'),
    MaxChi2 = cms.double(9.0),
    MaxDisplacement = cms.double(0.5),
    MaxSagitta = cms.double(2),
    MinPtForHitRecoveryInGluedDet = cms.double(1000000000000),
    MinimalTolerance = cms.double(0.5),
    appendToDataLabel = cms.string(''),
    clusterChargeCut = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutTight')
    ),
    nSigma = cms.double(3),
    pTChargeCutThreshold = cms.double(-1)
)


process.detachedQuadStepTrackCandidatesMkFitConfig = cms.ESProducer("MkFitIterationConfigESProducer",
    ComponentName = cms.string('detachedQuadStepTrackCandidatesMkFitConfig'),
    appendToDataLabel = cms.string(''),
    config = cms.FileInPath('RecoTracker/MkFit/data/mkfit-phase1-detachedQuadStep.json'),
    maxClusterSize = cms.uint32(8),
    minPt = cms.double(0)
)


process.detachedQuadStepTrajectoryCleanerBySharedHits = cms.ESProducer("TrajectoryCleanerESProducer",
    ComponentName = cms.string('detachedQuadStepTrajectoryCleanerBySharedHits'),
    ComponentType = cms.string('TrajectoryCleanerBySharedHits'),
    MissingHitPenalty = cms.double(20.0),
    ValidHitBonus = cms.double(5.0),
    allowSharedFirstHit = cms.bool(True),
    fractionShared = cms.double(0.13)
)


process.detachedTripletStepChi2Est = cms.ESProducer("Chi2ChargeMeasurementEstimatorESProducer",
    ComponentName = cms.string('detachedTripletStepChi2Est'),
    MaxChi2 = cms.double(9.0),
    MaxDisplacement = cms.double(0.5),
    MaxSagitta = cms.double(2),
    MinPtForHitRecoveryInGluedDet = cms.double(1000000000000),
    MinimalTolerance = cms.double(0.5),
    appendToDataLabel = cms.string(''),
    clusterChargeCut = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutTight')
    ),
    nSigma = cms.double(3),
    pTChargeCutThreshold = cms.double(-1)
)


process.detachedTripletStepTrackCandidatesMkFitConfig = cms.ESProducer("MkFitIterationConfigESProducer",
    ComponentName = cms.string('detachedTripletStepTrackCandidatesMkFitConfig'),
    appendToDataLabel = cms.string(''),
    config = cms.FileInPath('RecoTracker/MkFit/data/mkfit-phase1-detachedTripletStep.json'),
    maxClusterSize = cms.uint32(8),
    minPt = cms.double(0)
)


process.detachedTripletStepTrajectoryCleanerBySharedHits = cms.ESProducer("TrajectoryCleanerESProducer",
    ComponentName = cms.string('detachedTripletStepTrajectoryCleanerBySharedHits'),
    ComponentType = cms.string('TrajectoryCleanerBySharedHits'),
    MissingHitPenalty = cms.double(20.0),
    ValidHitBonus = cms.double(5.0),
    allowSharedFirstHit = cms.bool(True),
    fractionShared = cms.double(0.13)
)


process.displacedGeneralStepChi2Est = cms.ESProducer("Chi2ChargeMeasurementEstimatorESProducer",
    ComponentName = cms.string('displacedGeneralStepChi2Est'),
    MaxChi2 = cms.double(10.0),
    MaxDisplacement = cms.double(0.5),
    MaxSagitta = cms.double(2),
    MinPtForHitRecoveryInGluedDet = cms.double(1000000000000),
    MinimalTolerance = cms.double(0.5),
    appendToDataLabel = cms.string(''),
    clusterChargeCut = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutTight')
    ),
    nSigma = cms.double(3),
    pTChargeCutThreshold = cms.double(-1)
)


process.displacedGeneralStepClusterShapeHitFilter = cms.ESProducer("ClusterShapeHitFilterESProducer",
    ComponentName = cms.string('displacedGeneralStepClusterShapeHitFilter'),
    PixelShapeFile = cms.string('RecoPixelVertexing/PixelLowPtUtilities/data/pixelShapePhase1_noL1.par'),
    PixelShapeFileL1 = cms.string('RecoPixelVertexing/PixelLowPtUtilities/data/pixelShapePhase1_loose.par'),
    appendToDataLabel = cms.string(''),
    clusterChargeCut = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutTight')
    ),
    doPixelShapeCut = cms.bool(True),
    doStripShapeCut = cms.bool(False),
    isPhase2 = cms.bool(False)
)


process.displacedGeneralStepFitterSmoother = cms.ESProducer("KFFittingSmootherESProducer",
    BreakTrajWith2ConsecutiveMissing = cms.bool(True),
    ComponentName = cms.string('displacedGeneralStepFitterSmoother'),
    EstimateCut = cms.double(30),
    Fitter = cms.string('displacedGeneralStepRKFitter'),
    HighEtaSwitch = cms.double(5),
    LogPixelProbabilityCut = cms.double(0),
    MaxFractionOutliers = cms.double(0.3),
    MaxNumberOfOutliers = cms.int32(3),
    MinDof = cms.int32(2),
    MinNumberOfHits = cms.int32(8),
    MinNumberOfHitsHighEta = cms.int32(5),
    NoInvalidHitsBeginEnd = cms.bool(True),
    NoOutliersBeginEnd = cms.bool(False),
    RejectTracks = cms.bool(True),
    Smoother = cms.string('displacedGeneralStepRKSmoother'),
    appendToDataLabel = cms.string('')
)


process.displacedGeneralStepRKTrajectoryFitter = cms.ESProducer("KFTrajectoryFitterESProducer",
    ComponentName = cms.string('displacedGeneralStepRKFitter'),
    Estimator = cms.string('Chi2'),
    Propagator = cms.string('RungeKuttaTrackerPropagator'),
    RecoGeometry = cms.string('GlobalDetLayerGeometry'),
    Updator = cms.string('KFUpdator'),
    appendToDataLabel = cms.string(''),
    minHits = cms.int32(8)
)


process.displacedGeneralStepRKTrajectorySmoother = cms.ESProducer("KFTrajectorySmootherESProducer",
    ComponentName = cms.string('displacedGeneralStepRKSmoother'),
    Estimator = cms.string('Chi2'),
    Propagator = cms.string('RungeKuttaTrackerPropagator'),
    RecoGeometry = cms.string('GlobalDetLayerGeometry'),
    Updator = cms.string('KFUpdator'),
    appendToDataLabel = cms.string(''),
    errorRescaling = cms.double(10.0),
    minHits = cms.int32(8)
)


process.displacedGeneralStepTrajectoryCleanerBySharedHits = cms.ESProducer("TrajectoryCleanerESProducer",
    ComponentName = cms.string('displacedGeneralStepTrajectoryCleanerBySharedHits'),
    ComponentType = cms.string('TrajectoryCleanerBySharedHits'),
    MissingHitPenalty = cms.double(20.0),
    ValidHitBonus = cms.double(5.0),
    allowSharedFirstHit = cms.bool(True),
    fractionShared = cms.double(0.25)
)


process.duplicateTrackCandidatesChi2Est = cms.ESProducer("Chi2MeasurementEstimatorESProducer",
    ComponentName = cms.string('duplicateTrackCandidatesChi2Est'),
    MaxChi2 = cms.double(100),
    MaxDisplacement = cms.double(0.5),
    MaxSagitta = cms.double(2),
    MinPtForHitRecoveryInGluedDet = cms.double(1000000000000),
    MinimalTolerance = cms.double(0.5),
    appendToDataLabel = cms.string(''),
    nSigma = cms.double(3)
)


process.ecalSeverityLevel = cms.ESProducer("EcalSeverityLevelESProducer",
    dbstatusMask = cms.PSet(
        kBad = cms.vstring(
            'kNonRespondingIsolated',
            'kDeadVFE',
            'kDeadFE',
            'kNoDataNoTP'
        ),
        kGood = cms.vstring('kOk'),
        kProblematic = cms.vstring(
            'kDAC',
            'kNoLaser',
            'kNoisy',
            'kNNoisy',
            'kNNNoisy',
            'kNNNNoisy',
            'kNNNNNoisy',
            'kFixedG6',
            'kFixedG1',
            'kFixedG0'
        ),
        kRecovered = cms.vstring(),
        kTime = cms.vstring(),
        kWeird = cms.vstring()
    ),
    flagMask = cms.PSet(
        kBad = cms.vstring(
            'kFaultyHardware',
            'kDead',
            'kKilled'
        ),
        kGood = cms.vstring('kGood'),
        kProblematic = cms.vstring(
            'kPoorReco',
            'kPoorCalib',
            'kNoisy',
            'kSaturated'
        ),
        kRecovered = cms.vstring(
            'kLeadingEdgeRecovered',
            'kTowerRecovered'
        ),
        kTime = cms.vstring('kOutOfTime'),
        kWeird = cms.vstring(
            'kWeird',
            'kDiWeird'
        )
    ),
    timeThresh = cms.double(2.0)
)


process.ecalSimulationParametersEB = cms.ESProducer("EcalSimParametersESModule",
    appendToDataLabel = cms.string(''),
    fromDD4hep = cms.bool(True),
    name = cms.string('EcalHitsEB')
)


process.ecalSimulationParametersEE = cms.ESProducer("EcalSimParametersESModule",
    appendToDataLabel = cms.string(''),
    fromDD4hep = cms.bool(True),
    name = cms.string('EcalHitsEE')
)


process.ecalSimulationParametersES = cms.ESProducer("EcalSimParametersESModule",
    appendToDataLabel = cms.string(''),
    fromDD4hep = cms.bool(True),
    name = cms.string('EcalHitsES')
)


process.fakeForIdealAlignment = cms.ESProducer("FakeAlignmentProducer",
    appendToDataLabel = cms.string('fakeForIdeal')
)


process.generalDisplacedFlexibleKFFittingSmoother = cms.ESProducer("FlexibleKFFittingSmootherESProducer",
    ComponentName = cms.string('generalDisplacedFlexibleKFFittingSmoother'),
    appendToDataLabel = cms.string(''),
    looperFitter = cms.string('LooperFittingSmoother'),
    standardFitter = cms.string('displacedGeneralStepFitterSmoother')
)


process.hcalDDDRecConstants = cms.ESProducer("HcalDDDRecConstantsESModule",
    appendToDataLabel = cms.string('')
)


process.hcalDDDSimConstants = cms.ESProducer("HcalDDDSimConstantsESModule",
    appendToDataLabel = cms.string('')
)


process.hcalSimulationConstants = cms.ESProducer("HcalSimulationConstantsESModule",
    appendToDataLabel = cms.string('')
)


process.hcalSimulationParameters = cms.ESProducer("HcalSimParametersESModule",
    appendToDataLabel = cms.string(''),
    fromDD4hep = cms.bool(True)
)


process.hcalTopologyIdeal = cms.ESProducer("HcalTopologyIdealEP",
    Exclude = cms.untracked.string(''),
    MergePosition = cms.untracked.bool(False),
    appendToDataLabel = cms.string('')
)


process.hcal_db_producer = cms.ESProducer("HcalDbProducer",
    dump = cms.untracked.vstring(''),
    file = cms.untracked.string('')
)


process.highPtTripletStepChi2Est = cms.ESProducer("Chi2ChargeMeasurementEstimatorESProducer",
    ComponentName = cms.string('highPtTripletStepChi2Est'),
    MaxChi2 = cms.double(30),
    MaxDisplacement = cms.double(0.5),
    MaxSagitta = cms.double(2),
    MinPtForHitRecoveryInGluedDet = cms.double(1000000000000),
    MinimalTolerance = cms.double(0.5),
    appendToDataLabel = cms.string(''),
    clusterChargeCut = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutLoose')
    ),
    nSigma = cms.double(3),
    pTChargeCutThreshold = cms.double(15.0)
)


process.highPtTripletStepTrackCandidatesMkFitConfig = cms.ESProducer("MkFitIterationConfigESProducer",
    ComponentName = cms.string('highPtTripletStepTrackCandidatesMkFitConfig'),
    appendToDataLabel = cms.string(''),
    config = cms.FileInPath('RecoTracker/MkFit/data/mkfit-phase1-highPtTripletStep.json'),
    maxClusterSize = cms.uint32(8),
    minPt = cms.double(0)
)


process.highPtTripletStepTrajectoryCleanerBySharedHits = cms.ESProducer("TrajectoryCleanerESProducer",
    ComponentName = cms.string('highPtTripletStepTrajectoryCleanerBySharedHits'),
    ComponentType = cms.string('TrajectoryCleanerBySharedHits'),
    MissingHitPenalty = cms.double(20.0),
    ValidHitBonus = cms.double(5.0),
    allowSharedFirstHit = cms.bool(True),
    fractionShared = cms.double(0.16)
)


process.hitCollectorForCosmicDCSeeds = cms.ESProducer("Chi2MeasurementEstimatorESProducer",
    ComponentName = cms.string('hitCollectorForCosmicDCSeeds'),
    MaxChi2 = cms.double(100.0),
    MaxDisplacement = cms.double(100),
    MaxSagitta = cms.double(-1.0),
    MinPtForHitRecoveryInGluedDet = cms.double(1000000000000),
    MinimalTolerance = cms.double(0.5),
    appendToDataLabel = cms.string(''),
    nSigma = cms.double(4.0)
)


process.hitCollectorForOutInMuonSeeds = cms.ESProducer("Chi2MeasurementEstimatorESProducer",
    ComponentName = cms.string('hitCollectorForOutInMuonSeeds'),
    MaxChi2 = cms.double(100.0),
    MaxDisplacement = cms.double(0.5),
    MaxSagitta = cms.double(2),
    MinPtForHitRecoveryInGluedDet = cms.double(1000000000000),
    MinimalTolerance = cms.double(0.5),
    appendToDataLabel = cms.string(''),
    nSigma = cms.double(4.0)
)


process.idealForDigiCSCGeometry = cms.ESProducer("CSCGeometryESModule",
    alignmentsLabel = cms.string('fakeForIdeal'),
    appendToDataLabel = cms.string('idealForDigi'),
    applyAlignment = cms.bool(False),
    debugV = cms.untracked.bool(False),
    fromDD4hep = cms.bool(False),
    fromDDD = cms.bool(False),
    useCentreTIOffsets = cms.bool(False),
    useGangedStripsInME1a = cms.bool(False),
    useOnlyWiresInME1a = cms.bool(False),
    useRealWireGeometry = cms.bool(True)
)


process.idealForDigiDTGeometry = cms.ESProducer("DTGeometryESModule",
    DDDetector = cms.ESInputTag("",""),
    alignmentsLabel = cms.string('fakeForIdeal'),
    appendToDataLabel = cms.string('idealForDigi'),
    applyAlignment = cms.bool(False),
    attribute = cms.string('MuStructure'),
    fromDD4hep = cms.bool(False),
    fromDDD = cms.bool(False),
    value = cms.string('MuonBarrelDT')
)


process.idealForDigiTrackerGeometry = cms.ESProducer("TrackerDigiGeometryESModule",
    alignmentsLabel = cms.string('fakeForIdeal'),
    appendToDataLabel = cms.string('idealForDigi'),
    applyAlignment = cms.bool(False),
    fromDDD = cms.bool(False)
)


process.initialStepChi2Est = cms.ESProducer("Chi2ChargeMeasurementEstimatorESProducer",
    ComponentName = cms.string('initialStepChi2Est'),
    MaxChi2 = cms.double(30),
    MaxDisplacement = cms.double(0.5),
    MaxSagitta = cms.double(2),
    MinPtForHitRecoveryInGluedDet = cms.double(1000000000000),
    MinimalTolerance = cms.double(0.5),
    appendToDataLabel = cms.string(''),
    clusterChargeCut = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutLoose')
    ),
    nSigma = cms.double(3),
    pTChargeCutThreshold = cms.double(15.0)
)


process.initialStepChi2EstPreSplitting = cms.ESProducer("Chi2ChargeMeasurementEstimatorESProducer",
    ComponentName = cms.string('initialStepChi2EstPreSplitting'),
    MaxChi2 = cms.double(16.0),
    MaxDisplacement = cms.double(0.5),
    MaxSagitta = cms.double(2),
    MinPtForHitRecoveryInGluedDet = cms.double(1000000000000),
    MinimalTolerance = cms.double(0.5),
    appendToDataLabel = cms.string(''),
    clusterChargeCut = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutLoose')
    ),
    nSigma = cms.double(3),
    pTChargeCutThreshold = cms.double(-1)
)


process.initialStepTrackCandidatesMkFitConfig = cms.ESProducer("MkFitIterationConfigESProducer",
    ComponentName = cms.string('initialStepTrackCandidatesMkFitConfig'),
    appendToDataLabel = cms.string(''),
    config = cms.FileInPath('RecoTracker/MkFit/data/mkfit-phase1-initialStep.json'),
    maxClusterSize = cms.uint32(8),
    minPt = cms.double(0)
)


process.initialStepTrackCandidatesMkFitConfigPreSplitting = cms.ESProducer("MkFitIterationConfigESProducer",
    ComponentName = cms.string('initialStepTrackCandidatesMkFitConfigPreSplitting'),
    appendToDataLabel = cms.string(''),
    config = cms.FileInPath('RecoTracker/MkFit/data/mkfit-phase1-initialStep.json'),
    maxClusterSize = cms.uint32(8),
    minPt = cms.double(0)
)


process.jetCoreRegionalStepChi2Est = cms.ESProducer("Chi2MeasurementEstimatorESProducer",
    ComponentName = cms.string('jetCoreRegionalStepChi2Est'),
    MaxChi2 = cms.double(30),
    MaxDisplacement = cms.double(0.5),
    MaxSagitta = cms.double(2),
    MinPtForHitRecoveryInGluedDet = cms.double(1000000000000),
    MinimalTolerance = cms.double(0.5),
    appendToDataLabel = cms.string(''),
    nSigma = cms.double(3)
)


process.jetCoreRegionalStepDeepCoreTrajectoryCleaner = cms.ESProducer("TrajectoryCleanerESProducer",
    ComponentName = cms.string('jetCoreRegionalStepDeepCoreTrajectoryCleaner'),
    ComponentType = cms.string('TrajectoryCleanerBySharedHits'),
    MissingHitPenalty = cms.double(20.0),
    ValidHitBonus = cms.double(5.0),
    allowSharedFirstHit = cms.bool(True),
    fractionShared = cms.double(0.45)
)


process.lowPtQuadStepChi2Est = cms.ESProducer("Chi2ChargeMeasurementEstimatorESProducer",
    ComponentName = cms.string('lowPtQuadStepChi2Est'),
    MaxChi2 = cms.double(9.0),
    MaxDisplacement = cms.double(0.5),
    MaxSagitta = cms.double(2),
    MinPtForHitRecoveryInGluedDet = cms.double(1000000000000),
    MinimalTolerance = cms.double(0.5),
    appendToDataLabel = cms.string(''),
    clusterChargeCut = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutTight')
    ),
    nSigma = cms.double(3),
    pTChargeCutThreshold = cms.double(-1)
)


process.lowPtQuadStepTrackCandidatesMkFitConfig = cms.ESProducer("MkFitIterationConfigESProducer",
    ComponentName = cms.string('lowPtQuadStepTrackCandidatesMkFitConfig'),
    appendToDataLabel = cms.string(''),
    config = cms.FileInPath('RecoTracker/MkFit/data/mkfit-phase1-lowPtQuadStep.json'),
    maxClusterSize = cms.uint32(8),
    minPt = cms.double(0)
)


process.lowPtQuadStepTrajectoryCleanerBySharedHits = cms.ESProducer("TrajectoryCleanerESProducer",
    ComponentName = cms.string('lowPtQuadStepTrajectoryCleanerBySharedHits'),
    ComponentType = cms.string('TrajectoryCleanerBySharedHits'),
    MissingHitPenalty = cms.double(20.0),
    ValidHitBonus = cms.double(5.0),
    allowSharedFirstHit = cms.bool(True),
    fractionShared = cms.double(0.16)
)


process.lowPtTripletStepChi2Est = cms.ESProducer("Chi2ChargeMeasurementEstimatorESProducer",
    ComponentName = cms.string('lowPtTripletStepChi2Est'),
    MaxChi2 = cms.double(9.0),
    MaxDisplacement = cms.double(0.5),
    MaxSagitta = cms.double(2),
    MinPtForHitRecoveryInGluedDet = cms.double(1000000000000),
    MinimalTolerance = cms.double(0.5),
    appendToDataLabel = cms.string(''),
    clusterChargeCut = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutTight')
    ),
    nSigma = cms.double(3),
    pTChargeCutThreshold = cms.double(-1)
)


process.lowPtTripletStepTrackCandidatesMkFitConfig = cms.ESProducer("MkFitIterationConfigESProducer",
    ComponentName = cms.string('lowPtTripletStepTrackCandidatesMkFitConfig'),
    appendToDataLabel = cms.string(''),
    config = cms.FileInPath('RecoTracker/MkFit/data/mkfit-phase1-lowPtTripletStep.json'),
    maxClusterSize = cms.uint32(8),
    minPt = cms.double(0)
)


process.lowPtTripletStepTrajectoryCleanerBySharedHits = cms.ESProducer("TrajectoryCleanerESProducer",
    ComponentName = cms.string('lowPtTripletStepTrajectoryCleanerBySharedHits'),
    ComponentType = cms.string('TrajectoryCleanerBySharedHits'),
    MissingHitPenalty = cms.double(20.0),
    ValidHitBonus = cms.double(5.0),
    allowSharedFirstHit = cms.bool(True),
    fractionShared = cms.double(0.16)
)


process.mixedTripletStepChi2Est = cms.ESProducer("Chi2ChargeMeasurementEstimatorESProducer",
    ComponentName = cms.string('mixedTripletStepChi2Est'),
    MaxChi2 = cms.double(16.0),
    MaxDisplacement = cms.double(0.5),
    MaxSagitta = cms.double(2),
    MinPtForHitRecoveryInGluedDet = cms.double(1000000000000),
    MinimalTolerance = cms.double(0.5),
    appendToDataLabel = cms.string(''),
    clusterChargeCut = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutTight')
    ),
    nSigma = cms.double(3),
    pTChargeCutThreshold = cms.double(-1)
)


process.mixedTripletStepClusterShapeHitFilter = cms.ESProducer("ClusterShapeHitFilterESProducer",
    ComponentName = cms.string('mixedTripletStepClusterShapeHitFilter'),
    PixelShapeFile = cms.string('RecoPixelVertexing/PixelLowPtUtilities/data/pixelShapePhase1_noL1.par'),
    PixelShapeFileL1 = cms.string('RecoPixelVertexing/PixelLowPtUtilities/data/pixelShapePhase1_loose.par'),
    appendToDataLabel = cms.string(''),
    clusterChargeCut = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutTight')
    ),
    doPixelShapeCut = cms.bool(True),
    doStripShapeCut = cms.bool(True),
    isPhase2 = cms.bool(False)
)


process.mixedTripletStepPropagator = cms.ESProducer("PropagatorWithMaterialESProducer",
    ComponentName = cms.string('mixedTripletStepPropagator'),
    Mass = cms.double(0.105),
    MaxDPhi = cms.double(1.6),
    PropagationDirection = cms.string('alongMomentum'),
    SimpleMagneticField = cms.string(''),
    ptMin = cms.double(0.1),
    useRungeKutta = cms.bool(False)
)


process.mixedTripletStepPropagatorOpposite = cms.ESProducer("PropagatorWithMaterialESProducer",
    ComponentName = cms.string('mixedTripletStepPropagatorOpposite'),
    Mass = cms.double(0.105),
    MaxDPhi = cms.double(1.6),
    PropagationDirection = cms.string('oppositeToMomentum'),
    SimpleMagneticField = cms.string(''),
    ptMin = cms.double(0.1),
    useRungeKutta = cms.bool(False)
)


process.mixedTripletStepTrackCandidatesMkFitConfig = cms.ESProducer("MkFitIterationConfigESProducer",
    ComponentName = cms.string('mixedTripletStepTrackCandidatesMkFitConfig'),
    appendToDataLabel = cms.string(''),
    config = cms.FileInPath('RecoTracker/MkFit/data/mkfit-phase1-mixedTripletStep.json'),
    maxClusterSize = cms.uint32(8),
    minPt = cms.double(0)
)


process.mixedTripletStepTrajectoryCleanerBySharedHits = cms.ESProducer("TrajectoryCleanerESProducer",
    ComponentName = cms.string('mixedTripletStepTrajectoryCleanerBySharedHits'),
    ComponentType = cms.string('TrajectoryCleanerBySharedHits'),
    MissingHitPenalty = cms.double(20.0),
    ValidHitBonus = cms.double(5.0),
    allowSharedFirstHit = cms.bool(True),
    fractionShared = cms.double(0.11)
)


process.mkFitGeometryESProducer = cms.ESProducer("MkFitGeometryESProducer",
    appendToDataLabel = cms.string('')
)


process.multipleScatteringParametrisationMakerESProducer = cms.ESProducer("MultipleScatteringParametrisationMakerESProducer",
    appendToDataLabel = cms.string('')
)


process.muonGeometryConstants = cms.ESProducer("MuonGeometryConstantsESModule",
    appendToDataLabel = cms.string(''),
    fromDD4hep = cms.bool(True)
)


process.muonOffsetESProducer = cms.ESProducer("MuonOffsetESProducer",
    appendToDataLabel = cms.string(''),
    fromDD4hep = cms.bool(True),
    names = cms.vstring(
        'MuonCommonNumbering',
        'MuonBarrel',
        'MuonEndcap',
        'MuonBarrelWheels',
        'MuonBarrelStation1',
        'MuonBarrelStation2',
        'MuonBarrelStation3',
        'MuonBarrelStation4',
        'MuonBarrelSuperLayer',
        'MuonBarrelLayer',
        'MuonBarrelWire',
        'MuonRpcPlane1I',
        'MuonRpcPlane1O',
        'MuonRpcPlane2I',
        'MuonRpcPlane2O',
        'MuonRpcPlane3S',
        'MuonRpcPlane4',
        'MuonRpcChamberLeft',
        'MuonRpcChamberMiddle',
        'MuonRpcChamberRight',
        'MuonRpcEndcap1',
        'MuonRpcEndcap2',
        'MuonRpcEndcap3',
        'MuonRpcEndcap4',
        'MuonRpcEndcapSector',
        'MuonRpcEndcapChamberB1',
        'MuonRpcEndcapChamberB2',
        'MuonRpcEndcapChamberB3',
        'MuonRpcEndcapChamberC1',
        'MuonRpcEndcapChamberC2',
        'MuonRpcEndcapChamberC3',
        'MuonRpcEndcapChamberE1',
        'MuonRpcEndcapChamberE2',
        'MuonRpcEndcapChamberE3',
        'MuonRpcEndcapChamberF1',
        'MuonRpcEndcapChamberF2',
        'MuonRpcEndcapChamberF3',
        'MuonRpcEndcapChamberG1',
        'MuonRpcEndcapChamberH1',
        'MuonEndcapStation1',
        'MuonEndcapStation2',
        'MuonEndcapStation3',
        'MuonEndcapStation4',
        'MuonEndcapSubrings',
        'MuonEndcapSectors',
        'MuonEndcapLayers',
        'MuonEndcapRing1',
        'MuonEndcapRing2',
        'MuonEndcapRing3',
        'MuonEndcapRingA',
        'MuonGEMEndcap',
        'MuonGEMEndcap2',
        'MuonGEMSector',
        'MuonGEMChamber'
    )
)


process.muonSeededFittingSmootherWithOutliersRejectionAndRK = cms.ESProducer("KFFittingSmootherESProducer",
    BreakTrajWith2ConsecutiveMissing = cms.bool(False),
    ComponentName = cms.string('muonSeededFittingSmootherWithOutliersRejectionAndRK'),
    EstimateCut = cms.double(50.0),
    Fitter = cms.string('RKFitter'),
    HighEtaSwitch = cms.double(5),
    LogPixelProbabilityCut = cms.double(0),
    MaxFractionOutliers = cms.double(0.3),
    MaxNumberOfOutliers = cms.int32(3),
    MinDof = cms.int32(2),
    MinNumberOfHits = cms.int32(3),
    MinNumberOfHitsHighEta = cms.int32(5),
    NoInvalidHitsBeginEnd = cms.bool(True),
    NoOutliersBeginEnd = cms.bool(False),
    RejectTracks = cms.bool(True),
    Smoother = cms.string('RKSmoother'),
    appendToDataLabel = cms.string('')
)


process.muonSeededMeasurementEstimatorForInOut = cms.ESProducer("Chi2MeasurementEstimatorESProducer",
    ComponentName = cms.string('muonSeededMeasurementEstimatorForInOut'),
    MaxChi2 = cms.double(80.0),
    MaxDisplacement = cms.double(0.5),
    MaxSagitta = cms.double(-1.0),
    MinPtForHitRecoveryInGluedDet = cms.double(1000000000000),
    MinimalTolerance = cms.double(0.5),
    appendToDataLabel = cms.string(''),
    nSigma = cms.double(4.0)
)


process.muonSeededMeasurementEstimatorForOutIn = cms.ESProducer("Chi2MeasurementEstimatorESProducer",
    ComponentName = cms.string('muonSeededMeasurementEstimatorForOutIn'),
    MaxChi2 = cms.double(30),
    MaxDisplacement = cms.double(0.5),
    MaxSagitta = cms.double(-1.0),
    MinPtForHitRecoveryInGluedDet = cms.double(1000000000000),
    MinimalTolerance = cms.double(0.5),
    appendToDataLabel = cms.string(''),
    nSigma = cms.double(3)
)


process.muonSeededTrajectoryCleanerBySharedHits = cms.ESProducer("TrajectoryCleanerESProducer",
    ComponentName = cms.string('muonSeededTrajectoryCleanerBySharedHits'),
    ComponentType = cms.string('TrajectoryCleanerBySharedHits'),
    MissingHitPenalty = cms.double(1.0),
    ValidHitBonus = cms.double(1000.0),
    allowSharedFirstHit = cms.bool(True),
    fractionShared = cms.double(0.1)
)


process.myTTRHBuilderWithoutAngle4MixedPairs = cms.ESProducer("TkTransientTrackingRecHitBuilderESProducer",
    ComponentName = cms.string('TTRHBuilderWithoutAngle4MixedPairs'),
    ComputeCoarseLocalPositionFromDisk = cms.bool(False),
    Matcher = cms.string('StandardMatcher'),
    Phase2StripCPE = cms.string(''),
    PixelCPE = cms.string('PixelCPEGeneric'),
    StripCPE = cms.string('Fake'),
    appendToDataLabel = cms.string('')
)


process.myTTRHBuilderWithoutAngle4MixedTriplets = cms.ESProducer("TkTransientTrackingRecHitBuilderESProducer",
    ComponentName = cms.string('TTRHBuilderWithoutAngle4MixedTriplets'),
    ComputeCoarseLocalPositionFromDisk = cms.bool(False),
    Matcher = cms.string('StandardMatcher'),
    Phase2StripCPE = cms.string(''),
    PixelCPE = cms.string('PixelCPEGeneric'),
    StripCPE = cms.string('Fake'),
    appendToDataLabel = cms.string('')
)


process.myTTRHBuilderWithoutAngle4PixelPairs = cms.ESProducer("TkTransientTrackingRecHitBuilderESProducer",
    ComponentName = cms.string('TTRHBuilderWithoutAngle4PixelPairs'),
    ComputeCoarseLocalPositionFromDisk = cms.bool(False),
    Matcher = cms.string('StandardMatcher'),
    Phase2StripCPE = cms.string(''),
    PixelCPE = cms.string('PixelCPEGeneric'),
    StripCPE = cms.string('Fake'),
    appendToDataLabel = cms.string('')
)


process.myTTRHBuilderWithoutAngle4PixelTriplets = cms.ESProducer("TkTransientTrackingRecHitBuilderESProducer",
    ComponentName = cms.string('TTRHBuilderWithoutAngle4PixelTriplets'),
    ComputeCoarseLocalPositionFromDisk = cms.bool(False),
    Matcher = cms.string('StandardMatcher'),
    Phase2StripCPE = cms.string(''),
    PixelCPE = cms.string('PixelCPEGeneric'),
    StripCPE = cms.string('Fake'),
    appendToDataLabel = cms.string('')
)


process.navigationSchoolESProducer = cms.ESProducer("NavigationSchoolESProducer",
    ComponentName = cms.string('SimpleNavigationSchool'),
    SimpleMagneticField = cms.string('')
)


process.phase2StripCPEESProducer = cms.ESProducer("Phase2StripCPEESProducer",
    ComponentType = cms.string('Phase2StripCPE'),
    parameters = cms.PSet(
        LorentzAngle_DB = cms.bool(True),
        TanLorentzAnglePerTesla = cms.double(0.07)
    )
)


process.phase2StripCPEGeometricESProducer = cms.ESProducer("Phase2StripCPEESProducer",
    ComponentType = cms.string('Phase2StripCPEGeometric'),
    parameters = cms.PSet(

    )
)


process.pixelCPEFastESProducerPhase1 = cms.ESProducer("PixelCPEFastESProducerPhase1",
    Alpha2Order = cms.bool(True),
    ClusterProbComputationFlag = cms.int32(0),
    ComponentName = cms.string('PixelCPEFast'),
    EdgeClusterErrorX = cms.double(50),
    EdgeClusterErrorY = cms.double(85),
    LoadTemplatesFromDB = cms.bool(True),
    MagneticFieldRecord = cms.ESInputTag("",""),
    TruncatePixelCharge = cms.bool(True),
    UseErrorsFromTemplates = cms.bool(True),
    appendToDataLabel = cms.string(''),
    doLorentzFromAlignment = cms.bool(False),
    lAOffset = cms.double(0),
    lAWidthBPix = cms.double(0),
    lAWidthFPix = cms.double(0),
    useLAFromDB = cms.bool(True),
    useLAWidthFromDB = cms.bool(True),
    xerr_barrel_l1 = cms.vdouble(0.00115, 0.0012, 0.00088),
    xerr_barrel_l1_def = cms.double(0.0103),
    xerr_barrel_ln = cms.vdouble(0.00115, 0.0012, 0.00088),
    xerr_barrel_ln_def = cms.double(0.0103),
    xerr_endcap = cms.vdouble(0.002, 0.002),
    xerr_endcap_def = cms.double(0.002),
    yerr_barrel_l1 = cms.vdouble(
        0.00375, 0.0023, 0.0025, 0.0025, 0.0023,
        0.0023, 0.0021, 0.0021, 0.0024
    ),
    yerr_barrel_l1_def = cms.double(0.0021),
    yerr_barrel_ln = cms.vdouble(
        0.00375, 0.0023, 0.0025, 0.0025, 0.0023,
        0.0023, 0.0021, 0.0021, 0.0024
    ),
    yerr_barrel_ln_def = cms.double(0.0021),
    yerr_endcap = cms.vdouble(0.0021),
    yerr_endcap_def = cms.double(0.00075)
)


process.pixelCPEFastESProducerPhase2 = cms.ESProducer("PixelCPEFastESProducerPhase2",
    Alpha2Order = cms.bool(True),
    ClusterProbComputationFlag = cms.int32(0),
    ComponentName = cms.string('PixelCPEFastPhase2'),
    EdgeClusterErrorX = cms.double(50),
    EdgeClusterErrorY = cms.double(85),
    LoadTemplatesFromDB = cms.bool(True),
    MagneticFieldRecord = cms.ESInputTag("",""),
    TruncatePixelCharge = cms.bool(True),
    UseErrorsFromTemplates = cms.bool(True),
    appendToDataLabel = cms.string(''),
    doLorentzFromAlignment = cms.bool(False),
    lAOffset = cms.double(0),
    lAWidthBPix = cms.double(0),
    lAWidthFPix = cms.double(0),
    useLAFromDB = cms.bool(True),
    useLAWidthFromDB = cms.bool(True),
    xerr_barrel_l1 = cms.vdouble(0.00115, 0.0012, 0.00088),
    xerr_barrel_l1_def = cms.double(0.0103),
    xerr_barrel_ln = cms.vdouble(0.00115, 0.0012, 0.00088),
    xerr_barrel_ln_def = cms.double(0.0103),
    xerr_endcap = cms.vdouble(0.002, 0.002),
    xerr_endcap_def = cms.double(0.002),
    yerr_barrel_l1 = cms.vdouble(
        0.00375, 0.0023, 0.0025, 0.0025, 0.0023,
        0.0023, 0.0021, 0.0021, 0.0024
    ),
    yerr_barrel_l1_def = cms.double(0.0021),
    yerr_barrel_ln = cms.vdouble(
        0.00375, 0.0023, 0.0025, 0.0025, 0.0023,
        0.0023, 0.0021, 0.0021, 0.0024
    ),
    yerr_barrel_ln_def = cms.double(0.0021),
    yerr_endcap = cms.vdouble(0.0021),
    yerr_endcap_def = cms.double(0.00075)
)


process.pixelLessStepChi2Est = cms.ESProducer("Chi2ChargeMeasurementEstimatorESProducer",
    ComponentName = cms.string('pixelLessStepChi2Est'),
    MaxChi2 = cms.double(16.0),
    MaxDisplacement = cms.double(0.5),
    MaxSagitta = cms.double(2),
    MinPtForHitRecoveryInGluedDet = cms.double(1000000000000),
    MinimalTolerance = cms.double(0.5),
    appendToDataLabel = cms.string(''),
    clusterChargeCut = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutTight')
    ),
    nSigma = cms.double(3),
    pTChargeCutThreshold = cms.double(-1)
)


process.pixelLessStepClusterShapeHitFilter = cms.ESProducer("ClusterShapeHitFilterESProducer",
    ComponentName = cms.string('pixelLessStepClusterShapeHitFilter'),
    PixelShapeFile = cms.string('RecoPixelVertexing/PixelLowPtUtilities/data/pixelShapePhase1_noL1.par'),
    PixelShapeFileL1 = cms.string('RecoPixelVertexing/PixelLowPtUtilities/data/pixelShapePhase1_loose.par'),
    appendToDataLabel = cms.string(''),
    clusterChargeCut = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutTight')
    ),
    doPixelShapeCut = cms.bool(True),
    doStripShapeCut = cms.bool(False),
    isPhase2 = cms.bool(False)
)


process.pixelLessStepTrackCandidatesMkFitConfig = cms.ESProducer("MkFitIterationConfigESProducer",
    ComponentName = cms.string('pixelLessStepTrackCandidatesMkFitConfig'),
    appendToDataLabel = cms.string(''),
    config = cms.FileInPath('RecoTracker/MkFit/data/mkfit-phase1-pixelLessStep.json'),
    maxClusterSize = cms.uint32(8),
    minPt = cms.double(0)
)


process.pixelLessStepTrajectoryCleanerBySharedHits = cms.ESProducer("TrajectoryCleanerESProducer",
    ComponentName = cms.string('pixelLessStepTrajectoryCleanerBySharedHits'),
    ComponentType = cms.string('TrajectoryCleanerBySharedHits'),
    MissingHitPenalty = cms.double(20.0),
    ValidHitBonus = cms.double(5.0),
    allowSharedFirstHit = cms.bool(True),
    fractionShared = cms.double(0.11)
)


process.pixelPairStepChi2Est = cms.ESProducer("Chi2ChargeMeasurementEstimatorESProducer",
    ComponentName = cms.string('pixelPairStepChi2Est'),
    MaxChi2 = cms.double(9.0),
    MaxDisplacement = cms.double(0.5),
    MaxSagitta = cms.double(2),
    MinPtForHitRecoveryInGluedDet = cms.double(1000000000000),
    MinimalTolerance = cms.double(0.5),
    appendToDataLabel = cms.string(''),
    clusterChargeCut = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutLoose')
    ),
    nSigma = cms.double(3),
    pTChargeCutThreshold = cms.double(15.0)
)


process.pixelPairStepTrackCandidatesMkFitConfig = cms.ESProducer("MkFitIterationConfigESProducer",
    ComponentName = cms.string('pixelPairStepTrackCandidatesMkFitConfig'),
    appendToDataLabel = cms.string(''),
    config = cms.FileInPath('RecoTracker/MkFit/data/mkfit-phase1-pixelPairStep.json'),
    maxClusterSize = cms.uint32(8),
    minPt = cms.double(0)
)


process.pixelPairStepTrajectoryCleanerBySharedHits = cms.ESProducer("TrajectoryCleanerESProducer",
    ComponentName = cms.string('pixelPairStepTrajectoryCleanerBySharedHits'),
    ComponentType = cms.string('TrajectoryCleanerBySharedHits'),
    MissingHitPenalty = cms.double(20.0),
    ValidHitBonus = cms.double(5.0),
    allowSharedFirstHit = cms.bool(True),
    fractionShared = cms.double(0.095)
)


process.siPhase2RecHitMatcher = cms.ESProducer("SiPhase2RecHitMatcherESProducer",
    Algorithm = cms.string('VectorHitBuilderAlgorithm'),
    BarrelCut = cms.vdouble(
        0, 0.05, 0.06, 0.08, 0.09,
        0.12, 0.2
    ),
    CPE = cms.ESInputTag("","Phase2StripCPE"),
    Clusters = cms.string('siPhase2Clusters'),
    ComponentName = cms.string('SiPhase2VectorHitMatcher'),
    EndcapCut = cms.vdouble(
        0, 0.1, 0.1, 0.1, 0.1,
        0.1
    ),
    Phase2CPE_name = cms.string('Phase2StripCPE'),
    appendToDataLabel = cms.string(''),
    maxVectorHits = cms.int32(999999999),
    maxVectorHitsInAStack = cms.int32(999),
    offlinestubs = cms.string('vectorHits')
)


process.siPixel2DTemplateDBObjectESProducer = cms.ESProducer("SiPixel2DTemplateDBObjectESProducer")


process.siPixelGainCalibrationForHLTGPU = cms.ESProducer("SiPixelGainCalibrationForHLTGPUESProducer",
    appendToDataLabel = cms.string('')
)


process.siPixelQualityESProducer = cms.ESProducer("SiPixelQualityESProducer",
    ListOfRecordToMerge = cms.VPSet(
        cms.PSet(
            record = cms.string('SiPixelQualityFromDbRcd'),
            tag = cms.string('')
        ),
        cms.PSet(
            record = cms.string('SiPixelDetVOffRcd'),
            tag = cms.string('')
        )
    ),
    appendToDataLabel = cms.string(''),
    siPixelQualityLabel = cms.string(''),
    siPixelQualityLabel_RawToDigi = cms.string('')
)


process.siPixelROCsStatusAndMappingWrapperESProducer = cms.ESProducer("SiPixelROCsStatusAndMappingWrapperESProducer",
    CablingMapLabel = cms.string(''),
    ComponentName = cms.string(''),
    UseQualityInfo = cms.bool(False),
    appendToDataLabel = cms.string('')
)


process.siPixelTemplateDBObjectESProducer = cms.ESProducer("SiPixelTemplateDBObjectESProducer")


process.siStripBackPlaneCorrectionDepESProducer = cms.ESProducer("SiStripBackPlaneCorrectionDepESProducer",
    BackPlaneCorrectionDeconvMode = cms.PSet(
        label = cms.untracked.string('deconvolution'),
        record = cms.string('SiStripBackPlaneCorrectionRcd')
    ),
    BackPlaneCorrectionPeakMode = cms.PSet(
        label = cms.untracked.string('peak'),
        record = cms.string('SiStripBackPlaneCorrectionRcd')
    ),
    LatencyRecord = cms.PSet(
        label = cms.untracked.string(''),
        record = cms.string('SiStripLatencyRcd')
    )
)


process.siStripGainESProducer = cms.ESProducer("SiStripGainESProducer",
    APVGain = cms.VPSet(
        cms.PSet(
            Label = cms.untracked.string(''),
            NormalizationFactor = cms.untracked.double(1.0),
            Record = cms.string('SiStripApvGainRcd')
        ),
        cms.PSet(
            Label = cms.untracked.string(''),
            NormalizationFactor = cms.untracked.double(1.0),
            Record = cms.string('SiStripApvGain2Rcd')
        )
    ),
    AutomaticNormalization = cms.bool(False),
    appendToDataLabel = cms.string(''),
    printDebug = cms.untracked.bool(False)
)


process.siStripLorentzAngleDepESProducer = cms.ESProducer("SiStripLorentzAngleDepESProducer",
    LatencyRecord = cms.PSet(
        label = cms.untracked.string(''),
        record = cms.string('SiStripLatencyRcd')
    ),
    LorentzAngleDeconvMode = cms.PSet(
        label = cms.untracked.string('deconvolution'),
        record = cms.string('SiStripLorentzAngleRcd')
    ),
    LorentzAnglePeakMode = cms.PSet(
        label = cms.untracked.string('peak'),
        record = cms.string('SiStripLorentzAngleRcd')
    )
)


process.siStripQualityESProducer = cms.ESProducer("SiStripQualityESProducer",
    ListOfRecordToMerge = cms.VPSet(
        cms.PSet(
            record = cms.string('SiStripDetVOffRcd'),
            tag = cms.string('')
        ),
        cms.PSet(
            record = cms.string('SiStripDetCablingRcd'),
            tag = cms.string('')
        ),
        cms.PSet(
            record = cms.string('RunInfoRcd'),
            tag = cms.string('')
        ),
        cms.PSet(
            record = cms.string('SiStripBadChannelRcd'),
            tag = cms.string('')
        ),
        cms.PSet(
            record = cms.string('SiStripBadFiberRcd'),
            tag = cms.string('')
        ),
        cms.PSet(
            record = cms.string('SiStripBadModuleRcd'),
            tag = cms.string('')
        ),
        cms.PSet(
            record = cms.string('SiStripBadStripRcd'),
            tag = cms.string('')
        )
    ),
    PrintDebugOutput = cms.bool(False),
    ReduceGranularity = cms.bool(False),
    ThresholdForReducedGranularity = cms.double(0.3),
    UseEmptyRunInfo = cms.bool(False),
    appendToDataLabel = cms.string('')
)


process.sistripconn = cms.ESProducer("SiStripConnectivity")


process.stripCPEESProducer = cms.ESProducer("StripCPEESProducer",
    ComponentName = cms.string('stripCPE'),
    ComponentType = cms.string('SimpleStripCPE'),
    parameters = cms.PSet(

    )
)


process.templates = cms.ESProducer("PixelCPETemplateRecoESProducer",
    Alpha2Order = cms.bool(True),
    ClusterProbComputationFlag = cms.int32(0),
    ComponentName = cms.string('PixelCPETemplateReco'),
    LoadTemplatesFromDB = cms.bool(True),
    UseClusterSplitter = cms.bool(False),
    appendToDataLabel = cms.string(''),
    barrelTemplateID = cms.int32(0),
    directoryWithTemplates = cms.int32(0),
    doLorentzFromAlignment = cms.bool(False),
    forwardTemplateID = cms.int32(0),
    lAOffset = cms.double(0),
    lAWidthBPix = cms.double(0),
    lAWidthFPix = cms.double(0),
    speed = cms.int32(-2),
    useLAFromDB = cms.bool(True),
    useLAWidthFromDB = cms.bool(True)
)


process.templates2 = cms.ESProducer("PixelCPEClusterRepairESProducer",
    Alpha2Order = cms.bool(True),
    ClusterProbComputationFlag = cms.int32(0),
    ComponentName = cms.string('PixelCPEClusterRepair'),
    LoadTemplatesFromDB = cms.bool(True),
    MaxSizeMismatchInY = cms.double(0.3),
    MinChargeRatio = cms.double(0.8),
    Recommend2D = cms.vstring(
        'PXB 2',
        'PXB 3',
        'PXB 4'
    ),
    RunDamagedClusters = cms.bool(False),
    UseClusterSplitter = cms.bool(False),
    appendToDataLabel = cms.string(''),
    barrelTemplateID = cms.int32(0),
    directoryWithTemplates = cms.int32(0),
    doLorentzFromAlignment = cms.bool(False),
    forwardTemplateID = cms.int32(0),
    lAOffset = cms.double(0),
    lAWidthBPix = cms.double(0),
    lAWidthFPix = cms.double(0),
    speed = cms.int32(-2),
    useLAFromDB = cms.bool(True),
    useLAWidthFromDB = cms.bool(True)
)


process.templates2_speed0 = cms.ESProducer("PixelCPEClusterRepairESProducer",
    Alpha2Order = cms.bool(True),
    ClusterProbComputationFlag = cms.int32(0),
    ComponentName = cms.string('PixelCPEClusterRepairWithoutProbQ'),
    LoadTemplatesFromDB = cms.bool(True),
    MaxSizeMismatchInY = cms.double(0.3),
    MinChargeRatio = cms.double(0.8),
    Recommend2D = cms.vstring(
        'PXB 2',
        'PXB 3',
        'PXB 4'
    ),
    RunDamagedClusters = cms.bool(False),
    UseClusterSplitter = cms.bool(False),
    appendToDataLabel = cms.string(''),
    barrelTemplateID = cms.int32(0),
    directoryWithTemplates = cms.int32(0),
    doLorentzFromAlignment = cms.bool(False),
    forwardTemplateID = cms.int32(0),
    lAOffset = cms.double(0),
    lAWidthBPix = cms.double(0),
    lAWidthFPix = cms.double(0),
    speed = cms.int32(0),
    useLAFromDB = cms.bool(True),
    useLAWidthFromDB = cms.bool(True)
)


process.tkTransientTrackingRecHitBuilderESProducer = cms.ESProducer("TkTransientTrackingRecHitBuilderESProducer",
    ComponentName = cms.string('Fake'),
    ComputeCoarseLocalPositionFromDisk = cms.bool(False),
    Matcher = cms.string('Fake'),
    Phase2StripCPE = cms.string(''),
    PixelCPE = cms.string('Fake'),
    StripCPE = cms.string('Fake'),
    appendToDataLabel = cms.string('')
)


process.tobTecFlexibleKFFittingSmoother = cms.ESProducer("FlexibleKFFittingSmootherESProducer",
    ComponentName = cms.string('tobTecFlexibleKFFittingSmoother'),
    appendToDataLabel = cms.string(''),
    looperFitter = cms.string('tobTecStepFitterSmootherForLoopers'),
    standardFitter = cms.string('tobTecStepFitterSmoother')
)


process.tobTecStepChi2Est = cms.ESProducer("Chi2ChargeMeasurementEstimatorESProducer",
    ComponentName = cms.string('tobTecStepChi2Est'),
    MaxChi2 = cms.double(16.0),
    MaxDisplacement = cms.double(0.5),
    MaxSagitta = cms.double(2),
    MinPtForHitRecoveryInGluedDet = cms.double(1000000000000),
    MinimalTolerance = cms.double(0.5),
    appendToDataLabel = cms.string(''),
    clusterChargeCut = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutTight')
    ),
    nSigma = cms.double(3),
    pTChargeCutThreshold = cms.double(-1)
)


process.tobTecStepClusterShapeHitFilter = cms.ESProducer("ClusterShapeHitFilterESProducer",
    ComponentName = cms.string('tobTecStepClusterShapeHitFilter'),
    PixelShapeFile = cms.string('RecoPixelVertexing/PixelLowPtUtilities/data/pixelShapePhase1_noL1.par'),
    PixelShapeFileL1 = cms.string('RecoPixelVertexing/PixelLowPtUtilities/data/pixelShapePhase1_loose.par'),
    appendToDataLabel = cms.string(''),
    clusterChargeCut = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutTight')
    ),
    doPixelShapeCut = cms.bool(True),
    doStripShapeCut = cms.bool(False),
    isPhase2 = cms.bool(False)
)


process.tobTecStepFitterSmoother = cms.ESProducer("KFFittingSmootherESProducer",
    BreakTrajWith2ConsecutiveMissing = cms.bool(True),
    ComponentName = cms.string('tobTecStepFitterSmoother'),
    EstimateCut = cms.double(30),
    Fitter = cms.string('tobTecStepRKFitter'),
    HighEtaSwitch = cms.double(5),
    LogPixelProbabilityCut = cms.double(0),
    MaxFractionOutliers = cms.double(0.3),
    MaxNumberOfOutliers = cms.int32(3),
    MinDof = cms.int32(2),
    MinNumberOfHits = cms.int32(7),
    MinNumberOfHitsHighEta = cms.int32(5),
    NoInvalidHitsBeginEnd = cms.bool(True),
    NoOutliersBeginEnd = cms.bool(False),
    RejectTracks = cms.bool(True),
    Smoother = cms.string('tobTecStepRKSmoother'),
    appendToDataLabel = cms.string('')
)


process.tobTecStepFitterSmootherForLoopers = cms.ESProducer("KFFittingSmootherESProducer",
    BreakTrajWith2ConsecutiveMissing = cms.bool(True),
    ComponentName = cms.string('tobTecStepFitterSmootherForLoopers'),
    EstimateCut = cms.double(30),
    Fitter = cms.string('tobTecStepRKFitterForLoopers'),
    HighEtaSwitch = cms.double(5),
    LogPixelProbabilityCut = cms.double(0),
    MaxFractionOutliers = cms.double(0.3),
    MaxNumberOfOutliers = cms.int32(3),
    MinDof = cms.int32(2),
    MinNumberOfHits = cms.int32(7),
    MinNumberOfHitsHighEta = cms.int32(5),
    NoInvalidHitsBeginEnd = cms.bool(True),
    NoOutliersBeginEnd = cms.bool(False),
    RejectTracks = cms.bool(True),
    Smoother = cms.string('tobTecStepRKSmootherForLoopers'),
    appendToDataLabel = cms.string('')
)


process.tobTecStepRKTrajectoryFitter = cms.ESProducer("KFTrajectoryFitterESProducer",
    ComponentName = cms.string('tobTecStepRKFitter'),
    Estimator = cms.string('Chi2'),
    Propagator = cms.string('RungeKuttaTrackerPropagator'),
    RecoGeometry = cms.string('GlobalDetLayerGeometry'),
    Updator = cms.string('KFUpdator'),
    appendToDataLabel = cms.string(''),
    minHits = cms.int32(7)
)


process.tobTecStepRKTrajectoryFitterForLoopers = cms.ESProducer("KFTrajectoryFitterESProducer",
    ComponentName = cms.string('tobTecStepRKFitterForLoopers'),
    Estimator = cms.string('Chi2'),
    Propagator = cms.string('PropagatorWithMaterialForLoopers'),
    RecoGeometry = cms.string('GlobalDetLayerGeometry'),
    Updator = cms.string('KFUpdator'),
    appendToDataLabel = cms.string(''),
    minHits = cms.int32(7)
)


process.tobTecStepRKTrajectorySmoother = cms.ESProducer("KFTrajectorySmootherESProducer",
    ComponentName = cms.string('tobTecStepRKSmoother'),
    Estimator = cms.string('Chi2'),
    Propagator = cms.string('RungeKuttaTrackerPropagator'),
    RecoGeometry = cms.string('GlobalDetLayerGeometry'),
    Updator = cms.string('KFUpdator'),
    appendToDataLabel = cms.string(''),
    errorRescaling = cms.double(10.0),
    minHits = cms.int32(7)
)


process.tobTecStepRKTrajectorySmootherForLoopers = cms.ESProducer("KFTrajectorySmootherESProducer",
    ComponentName = cms.string('tobTecStepRKSmootherForLoopers'),
    Estimator = cms.string('Chi2'),
    Propagator = cms.string('PropagatorWithMaterialForLoopers'),
    RecoGeometry = cms.string('GlobalDetLayerGeometry'),
    Updator = cms.string('KFUpdator'),
    appendToDataLabel = cms.string(''),
    errorRescaling = cms.double(10.0),
    minHits = cms.int32(7)
)


process.tobTecStepTrackCandidatesMkFitConfig = cms.ESProducer("MkFitIterationConfigESProducer",
    ComponentName = cms.string('tobTecStepTrackCandidatesMkFitConfig'),
    appendToDataLabel = cms.string(''),
    config = cms.FileInPath('RecoTracker/MkFit/data/mkfit-phase1-tobTecStep.json'),
    maxClusterSize = cms.uint32(8),
    minPt = cms.double(0)
)


process.tobTecStepTrajectoryCleanerBySharedHits = cms.ESProducer("TrajectoryCleanerESProducer",
    ComponentName = cms.string('tobTecStepTrajectoryCleanerBySharedHits'),
    ComponentType = cms.string('TrajectoryCleanerBySharedHits'),
    MissingHitPenalty = cms.double(20.0),
    ValidHitBonus = cms.double(5.0),
    allowSharedFirstHit = cms.bool(True),
    fractionShared = cms.double(0.09)
)


process.trackAlgoPriorityOrder = cms.ESProducer("TrackAlgoPriorityOrderESProducer",
    ComponentName = cms.string('trackAlgoPriorityOrder'),
    algoOrder = cms.vstring(
        'initialStep',
        'lowPtQuadStep',
        'highPtTripletStep',
        'lowPtTripletStep',
        'detachedQuadStep',
        'detachedTripletStep',
        'pixelPairStep',
        'mixedTripletStep',
        'pixelLessStep',
        'tobTecStep',
        'jetCoreRegionalStep',
        'muonSeededStepInOut',
        'muonSeededStepOutIn'
    ),
    appendToDataLabel = cms.string('')
)


process.trackSelectionTf = cms.ESProducer("TfGraphDefProducer",
    ComponentName = cms.string('trackSelectionTf'),
    FileName = cms.FileInPath('RecoTracker/FinalTrackSelectors/data/TrackTfClassifier/MkFit_Run3_12_5_0_pre5.pb'),
    appendToDataLabel = cms.string('')
)


process.trackSelectionTfPLess = cms.ESProducer("TfGraphDefProducer",
    ComponentName = cms.string('trackSelectionTfPLess'),
    FileName = cms.FileInPath('RecoTracker/FinalTrackSelectors/data/TrackTfClassifier/MkFitPixelLessOnly_Run3_12_5_0_pre5.pb'),
    appendToDataLabel = cms.string('')
)


process.trackSelectionTf_CKF = cms.ESProducer("TfGraphDefProducer",
    ComponentName = cms.string('trackSelectionTf_CKF'),
    FileName = cms.FileInPath('RecoTracker/FinalTrackSelectors/data/TrackTfClassifier/CKF_Run3_12_5_0_pre5.pb'),
    appendToDataLabel = cms.string('')
)


process.trackerGeometryDB = cms.ESProducer("TrackerDigiGeometryESModule",
    alignmentsLabel = cms.string(''),
    appendToDataLabel = cms.string(''),
    applyAlignment = cms.bool(True),
    fromDDD = cms.bool(False)
)


process.trackerNumberingGeometryDB = cms.ESProducer("TrackerGeometricDetESModule",
    appendToDataLabel = cms.string(''),
    fromDD4hep = cms.bool(False),
    fromDDD = cms.bool(False)
)


process.trackerTopology = cms.ESProducer("TrackerTopologyEP",
    appendToDataLabel = cms.string('')
)


process.trajectoryCleanerBySharedHits = cms.ESProducer("TrajectoryCleanerESProducer",
    ComponentName = cms.string('TrajectoryCleanerBySharedHits'),
    ComponentType = cms.string('TrajectoryCleanerBySharedHits'),
    MissingHitPenalty = cms.double(20.0),
    ValidHitBonus = cms.double(5.0),
    allowSharedFirstHit = cms.bool(True),
    fractionShared = cms.double(0.19)
)


process.ttrhbwor = cms.ESProducer("TkTransientTrackingRecHitBuilderESProducer",
    ComponentName = cms.string('WithoutRefit'),
    ComputeCoarseLocalPositionFromDisk = cms.bool(False),
    Matcher = cms.string('Fake'),
    Phase2StripCPE = cms.string(''),
    PixelCPE = cms.string('Fake'),
    StripCPE = cms.string('Fake'),
    appendToDataLabel = cms.string('')
)


process.ttrhbwr = cms.ESProducer("TkTransientTrackingRecHitBuilderESProducer",
    ComponentName = cms.string('WithTrackAngle'),
    ComputeCoarseLocalPositionFromDisk = cms.bool(False),
    Matcher = cms.string('StandardMatcher'),
    Phase2StripCPE = cms.string(''),
    PixelCPE = cms.string('PixelCPEGeneric'),
    StripCPE = cms.string('StripCPEfromTrackAngle'),
    appendToDataLabel = cms.string('')
)


process.DDDetectorESProducerFromDB = cms.ESSource("DDDetectorESProducer",
    appendToDataLabel = cms.string(''),
    fromDB = cms.bool(True),
    label = cms.string('Extended'),
    rootDDName = cms.string('cms:OCMS')
)


process.GlobalTag = cms.ESSource("PoolDBESSource",
    DBParameters = cms.PSet(
        authenticationPath = cms.untracked.string(''),
        authenticationSystem = cms.untracked.int32(0),
        connectionTimeout = cms.untracked.int32(0),
        messageLevel = cms.untracked.int32(0),
        security = cms.untracked.string('')
    ),
    DumpStat = cms.untracked.bool(False),
    ReconnectEachRun = cms.untracked.bool(False),
    RefreshAlways = cms.untracked.bool(False),
    RefreshEachRun = cms.untracked.bool(False),
    RefreshOpenIOVs = cms.untracked.bool(False),
    connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
    globaltag = cms.string('130X_mcRun3_2023_realistic_postBPix_v5'),
    pfnPostfix = cms.untracked.string(''),
    pfnPrefix = cms.untracked.string(''),
    snapshotTime = cms.string(''),
    toGet = cms.VPSet()
)


process.HcalTimeSlewEP = cms.ESSource("HcalTimeSlewEP",
    appendToDataLabel = cms.string('HBHE'),
    timeSlewParametersM2 = cms.VPSet(
        cms.PSet(
            slope = cms.double(-3.178648),
            tmax = cms.double(16.0),
            tzero = cms.double(23.960177)
        ),
        cms.PSet(
            slope = cms.double(-1.5610227),
            tmax = cms.double(10.0),
            tzero = cms.double(11.977461)
        ),
        cms.PSet(
            slope = cms.double(-1.075824),
            tmax = cms.double(6.25),
            tzero = cms.double(9.109694)
        )
    ),
    timeSlewParametersM3 = cms.VPSet(
        cms.PSet(
            cap = cms.double(6.0),
            tspar0 = cms.double(12.2999),
            tspar0_siPM = cms.double(0.0),
            tspar1 = cms.double(-2.19142),
            tspar1_siPM = cms.double(0.0),
            tspar2 = cms.double(0.0),
            tspar2_siPM = cms.double(0.0)
        ),
        cms.PSet(
            cap = cms.double(6.0),
            tspar0 = cms.double(15.5),
            tspar0_siPM = cms.double(0.0),
            tspar1 = cms.double(-3.2),
            tspar1_siPM = cms.double(0.0),
            tspar2 = cms.double(32.0),
            tspar2_siPM = cms.double(0.0)
        ),
        cms.PSet(
            cap = cms.double(6.0),
            tspar0 = cms.double(12.2999),
            tspar0_siPM = cms.double(0.0),
            tspar1 = cms.double(-2.19142),
            tspar1_siPM = cms.double(0.0),
            tspar2 = cms.double(0.0),
            tspar2_siPM = cms.double(0.0)
        ),
        cms.PSet(
            cap = cms.double(6.0),
            tspar0 = cms.double(12.2999),
            tspar0_siPM = cms.double(0.0),
            tspar1 = cms.double(-2.19142),
            tspar1_siPM = cms.double(0.0),
            tspar2 = cms.double(0.0),
            tspar2_siPM = cms.double(0.0)
        )
    )
)


process.eegeom = cms.ESSource("EmptyESSource",
    firstValid = cms.vuint32(1),
    iovIsRunNotTime = cms.bool(True),
    recordName = cms.string('EcalMappingRcd')
)


process.es_hardcode = cms.ESSource("HcalHardcodeCalibrations",
    GainWidthsForTrigPrims = cms.bool(False),
    HBRecalibration = cms.bool(False),
    HBmeanenergies = cms.FileInPath('CalibCalorimetry/HcalPlugins/data/meanenergiesHB.txt'),
    HBreCalibCutoff = cms.double(100.0),
    HERecalibration = cms.bool(False),
    HEmeanenergies = cms.FileInPath('CalibCalorimetry/HcalPlugins/data/meanenergiesHE.txt'),
    HEreCalibCutoff = cms.double(100.0),
    HFRecalParameterBlock = cms.PSet(
        HFdepthOneParameterA = cms.vdouble(
            0.004123, 0.00602, 0.008201, 0.010489, 0.013379,
            0.016997, 0.021464, 0.027371, 0.034195, 0.044807,
            0.058939, 0.125497
        ),
        HFdepthOneParameterB = cms.vdouble(
            -4e-06, -2e-06, 0.0, 4e-06, 1.5e-05,
            2.6e-05, 6.3e-05, 8.4e-05, 0.00016, 0.000107,
            0.000425, 0.000209
        ),
        HFdepthTwoParameterA = cms.vdouble(
            0.002861, 0.004168, 0.0064, 0.008388, 0.011601,
            0.014425, 0.018633, 0.023232, 0.028274, 0.035447,
            0.051579, 0.086593
        ),
        HFdepthTwoParameterB = cms.vdouble(
            -2e-06, -0.0, -7e-06, -6e-06, -2e-06,
            1e-06, 1.9e-05, 3.1e-05, 6.7e-05, 1.2e-05,
            0.000157, -3e-06
        )
    ),
    HFRecalibration = cms.bool(False),
    SiPMCharacteristics = cms.VPSet(
        cms.PSet(
            crosstalk = cms.double(0.0),
            nonlin1 = cms.double(1.0),
            nonlin2 = cms.double(0.0),
            nonlin3 = cms.double(0.0),
            pixels = cms.int32(36000)
        ),
        cms.PSet(
            crosstalk = cms.double(0.0),
            nonlin1 = cms.double(1.0),
            nonlin2 = cms.double(0.0),
            nonlin3 = cms.double(0.0),
            pixels = cms.int32(2500)
        ),
        cms.PSet(
            crosstalk = cms.double(0.17),
            nonlin1 = cms.double(1.00985),
            nonlin2 = cms.double(7.84089e-06),
            nonlin3 = cms.double(2.86282e-10),
            pixels = cms.int32(27370)
        ),
        cms.PSet(
            crosstalk = cms.double(0.196),
            nonlin1 = cms.double(1.00546),
            nonlin2 = cms.double(6.40239e-06),
            nonlin3 = cms.double(1.27011e-10),
            pixels = cms.int32(38018)
        ),
        cms.PSet(
            crosstalk = cms.double(0.17),
            nonlin1 = cms.double(1.00985),
            nonlin2 = cms.double(7.84089e-06),
            nonlin3 = cms.double(2.86282e-10),
            pixels = cms.int32(27370)
        ),
        cms.PSet(
            crosstalk = cms.double(0.196),
            nonlin1 = cms.double(1.00546),
            nonlin2 = cms.double(6.40239e-06),
            nonlin3 = cms.double(1.27011e-10),
            pixels = cms.int32(38018)
        ),
        cms.PSet(
            crosstalk = cms.double(0.0),
            nonlin1 = cms.double(1.0),
            nonlin2 = cms.double(0.0),
            nonlin3 = cms.double(0.0),
            pixels = cms.int32(0)
        )
    ),
    hb = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.19),
        gainWidth = cms.vdouble(0.0),
        mcShape = cms.int32(125),
        noiseCorrelation = cms.vdouble(0.0),
        pedestal = cms.double(3.285),
        pedestalWidth = cms.double(0.809),
        photoelectronsToAnalog = cms.double(0.3305),
        qieOffset = cms.vdouble(-0.49, 1.8, 7.2, 37.9),
        qieSlope = cms.vdouble(0.912, 0.917, 0.922, 0.923),
        qieType = cms.int32(0),
        recoShape = cms.int32(105),
        zsThreshold = cms.int32(8)
    ),
    hbUpgrade = cms.PSet(
        darkCurrent = cms.vdouble(0.01, 0.015),
        doRadiationDamage = cms.bool(True),
        gain = cms.vdouble(0.0006252),
        gainWidth = cms.vdouble(0),
        mcShape = cms.int32(206),
        noiseCorrelation = cms.vdouble(0.26, 0.254),
        pedestal = cms.double(17.3),
        pedestalWidth = cms.double(1.5),
        photoelectronsToAnalog = cms.double(40.0),
        qieOffset = cms.vdouble(0.0, 0.0, 0.0, 0.0),
        qieSlope = cms.vdouble(0.05376, 0.05376, 0.05376, 0.05376),
        qieType = cms.int32(2),
        radiationDamage = cms.PSet(
            depVsNeutrons = cms.vdouble(5.543e-10, 8.012e-10),
            depVsTemp = cms.double(0.0631),
            intlumiOffset = cms.double(150),
            intlumiToNeutrons = cms.double(367000000.0),
            temperatureBase = cms.double(20),
            temperatureNew = cms.double(-5)
        ),
        recoShape = cms.int32(208),
        zsThreshold = cms.int32(16)
    ),
    he = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.23),
        gainWidth = cms.vdouble(0),
        mcShape = cms.int32(125),
        noiseCorrelation = cms.vdouble(0.0),
        pedestal = cms.double(3.163),
        pedestalWidth = cms.double(0.9698),
        photoelectronsToAnalog = cms.double(0.3305),
        qieOffset = cms.vdouble(-0.38, 2.0, 7.6, 39.6),
        qieSlope = cms.vdouble(0.912, 0.916, 0.92, 0.922),
        qieType = cms.int32(0),
        recoShape = cms.int32(105),
        zsThreshold = cms.int32(9)
    ),
    heUpgrade = cms.PSet(
        darkCurrent = cms.vdouble(0.01, 0.015),
        doRadiationDamage = cms.bool(True),
        gain = cms.vdouble(0.0006252),
        gainWidth = cms.vdouble(0),
        mcShape = cms.int32(206),
        noiseCorrelation = cms.vdouble(0.26, 0.254),
        pedestal = cms.double(17.3),
        pedestalWidth = cms.double(1.5),
        photoelectronsToAnalog = cms.double(40.0),
        qieOffset = cms.vdouble(0.0, 0.0, 0.0, 0.0),
        qieSlope = cms.vdouble(0.05376, 0.05376, 0.05376, 0.05376),
        qieType = cms.int32(2),
        radiationDamage = cms.PSet(
            depVsNeutrons = cms.vdouble(5.543e-10, 8.012e-10),
            depVsTemp = cms.double(0.0631),
            intlumiOffset = cms.double(75),
            intlumiToNeutrons = cms.double(29200000.0),
            temperatureBase = cms.double(20),
            temperatureNew = cms.double(5)
        ),
        recoShape = cms.int32(208),
        zsThreshold = cms.int32(16)
    ),
    hf = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.14, 0.135),
        gainWidth = cms.vdouble(0.0, 0.0),
        mcShape = cms.int32(301),
        noiseCorrelation = cms.vdouble(0.0),
        pedestal = cms.double(9.354),
        pedestalWidth = cms.double(2.516),
        photoelectronsToAnalog = cms.double(0.0),
        qieOffset = cms.vdouble(-0.87, 1.4, 7.8, -29.6),
        qieSlope = cms.vdouble(0.359, 0.358, 0.36, 0.367),
        qieType = cms.int32(0),
        recoShape = cms.int32(301),
        zsThreshold = cms.int32(-9999)
    ),
    hfUpgrade = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.14, 0.135),
        gainWidth = cms.vdouble(0.0, 0.0),
        mcShape = cms.int32(301),
        noiseCorrelation = cms.vdouble(0.0),
        pedestal = cms.double(13.33),
        pedestalWidth = cms.double(3.33),
        photoelectronsToAnalog = cms.double(0.0),
        qieOffset = cms.vdouble(0.0697, -0.7405, 12.38, -671.9),
        qieSlope = cms.vdouble(0.297, 0.298, 0.298, 0.313),
        qieType = cms.int32(1),
        recoShape = cms.int32(301),
        zsThreshold = cms.int32(-9999)
    ),
    ho = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.006, 0.0087),
        gainWidth = cms.vdouble(0.0, 0.0),
        mcShape = cms.int32(201),
        noiseCorrelation = cms.vdouble(0.0),
        pedestal = cms.double(12.06),
        pedestalWidth = cms.double(0.6285),
        photoelectronsToAnalog = cms.double(4.0),
        qieOffset = cms.vdouble(-0.44, 1.4, 7.1, 38.5),
        qieSlope = cms.vdouble(0.907, 0.915, 0.92, 0.921),
        qieType = cms.int32(0),
        recoShape = cms.int32(201),
        zsThreshold = cms.int32(24)
    ),
    iLumi = cms.double(-1.0),
    killHE = cms.bool(False),
    testHEPlan1 = cms.bool(False),
    testHFQIE10 = cms.bool(False),
    toGet = cms.untracked.vstring('GainWidths'),
    useHBUpgrade = cms.bool(True),
    useHEUpgrade = cms.bool(True),
    useHFUpgrade = cms.bool(True),
    useHOUpgrade = cms.bool(True),
    useIeta18depth1 = cms.bool(False),
    useLayer0Weight = cms.bool(True)
)


process.essourceEcalSev = cms.ESSource("EmptyESSource",
    firstValid = cms.vuint32(1),
    iovIsRunNotTime = cms.bool(True),
    recordName = cms.string('EcalSeverityLevelAlgoRcd')
)


process.trackdnn_source = cms.ESSource("EmptyESSource",
    firstValid = cms.vuint32(1),
    iovIsRunNotTime = cms.bool(True),
    recordName = cms.string('TfGraphRecord')
)


process.prefer("GlobalTag")

process.prefer("es_hardcode")

process.Conv2StepTask = cms.Task(process.conv2Clusters, process.conv2LayerPairs, process.conv2StepSelector, process.conv2StepTracks, process.conv2TrackCandidates, process.photonConvTrajSeedFromQuadruplets)


process.ConvStepTask = cms.Task(process.convClusters, process.convLayerPairs, process.convStepSelector, process.convStepTracks, process.convTrackCandidates, process.photonConvTrajSeedFromSingleLeg)


process.DetachedQuadStepTask = cms.Task(process.detachedQuadStep, process.detachedQuadStepClusters, process.detachedQuadStepHitDoublets, process.detachedQuadStepHitQuadruplets, process.detachedQuadStepSeedLayers, process.detachedQuadStepSeeds, process.detachedQuadStepTrackCandidates, process.detachedQuadStepTrackCandidatesMkFit, process.detachedQuadStepTrackCandidatesMkFitConfig, process.detachedQuadStepTrackCandidatesMkFitSeeds, process.detachedQuadStepTrackingRegions, process.detachedQuadStepTracks)


process.DetachedTripletStepTask = cms.Task(process.detachedTripletStep, process.detachedTripletStepClassifier1, process.detachedTripletStepClassifier2, process.detachedTripletStepClusters, process.detachedTripletStepHitDoublets, process.detachedTripletStepHitTriplets, process.detachedTripletStepSeedLayers, process.detachedTripletStepSeeds, process.detachedTripletStepTrackCandidates, process.detachedTripletStepTrackCandidatesMkFit, process.detachedTripletStepTrackCandidatesMkFitConfig, process.detachedTripletStepTrackCandidatesMkFitSeeds, process.detachedTripletStepTrackingRegions, process.detachedTripletStepTracks)


process.DisplacedGeneralStepTask = cms.Task(process.displacedGeneralStep, process.displacedGeneralStepClassifier1, process.displacedGeneralStepClassifier2, process.displacedGeneralStepClusters, process.displacedGeneralStepHitDoublets, process.displacedGeneralStepHitTriplets, process.displacedGeneralStepSeedLayers, process.displacedGeneralStepSeeds, process.displacedGeneralStepTrackCandidates, process.displacedGeneralStepTrackingRegions, process.displacedGeneralStepTracks)


process.HighPtTripletStepTask = cms.Task(process.highPtTripletStep, process.highPtTripletStepClusters, process.highPtTripletStepHitDoublets, process.highPtTripletStepHitTriplets, process.highPtTripletStepSeedLayers, process.highPtTripletStepSeeds, process.highPtTripletStepTrackCandidates, process.highPtTripletStepTrackCandidatesMkFit, process.highPtTripletStepTrackCandidatesMkFitConfig, process.highPtTripletStepTrackCandidatesMkFitSeeds, process.highPtTripletStepTrackingRegions, process.highPtTripletStepTracks)


process.InitialStepPreSplittingTask = cms.Task(process.MeasurementTrackerEvent, process.ak4CaloJetsForTrkPreSplitting, process.caloTowerForTrkPreSplitting, process.firstStepPrimaryVerticesPreSplitting, process.initialStepHitDoubletsPreSplitting, process.initialStepHitQuadrupletsPreSplitting, process.initialStepSeedLayersPreSplitting, process.initialStepSeedsPreSplitting, process.initialStepTrackCandidatesMkFitConfigPreSplitting, process.initialStepTrackCandidatesMkFitPreSplitting, process.initialStepTrackCandidatesMkFitSeedsPreSplitting, process.initialStepTrackCandidatesPreSplitting, process.initialStepTrackRefsForJetsPreSplitting, process.initialStepTrackingRegionsPreSplitting, process.initialStepTracksPreSplitting, process.jetsForCoreTrackingPreSplitting, process.mkFitEventOfHitsPreSplitting, process.mkFitGeometryESProducer, process.mkFitSiPixelHitsPreSplitting, process.mkFitSiStripHits, process.siPixelClusterShapeCache, process.siPixelClusters, process.siPixelRecHits, process.trackerClusterCheckPreSplitting)


process.JetCoreRegionalStepBarrelTask = cms.Task(process.firstStepGoodPrimaryVertices, process.jetCoreRegionalStepBarrel, process.jetCoreRegionalStepBarrelTrackCandidates, process.jetCoreRegionalStepBarrelTracks, process.jetCoreRegionalStepSeedLayers, process.jetCoreRegionalStepSeedsBarrel, process.jetsForCoreTrackingBarrel)


process.JetCoreRegionalStepEndcapTask = cms.Task(process.firstStepGoodPrimaryVertices, process.jetCoreRegionalStepEndcap, process.jetCoreRegionalStepEndcapHitDoublets, process.jetCoreRegionalStepEndcapTrackCandidates, process.jetCoreRegionalStepEndcapTrackingRegions, process.jetCoreRegionalStepEndcapTracks, process.jetCoreRegionalStepSeedLayers, process.jetCoreRegionalStepSeedsEndcap, process.jetsForCoreTrackingEndcap)


process.JetCoreRegionalStepTask = cms.Task(process.firstStepGoodPrimaryVertices, process.jetCoreRegionalStep, process.jetCoreRegionalStepHitDoublets, process.jetCoreRegionalStepSeedLayers, process.jetCoreRegionalStepSeeds, process.jetCoreRegionalStepTrackCandidates, process.jetCoreRegionalStepTrackingRegions, process.jetCoreRegionalStepTracks, process.jetsForCoreTracking)


process.LowPtQuadStepTask = cms.Task(process.lowPtQuadStep, process.lowPtQuadStepClusters, process.lowPtQuadStepHitDoublets, process.lowPtQuadStepHitQuadruplets, process.lowPtQuadStepSeedLayers, process.lowPtQuadStepSeeds, process.lowPtQuadStepTrackCandidates, process.lowPtQuadStepTrackingRegions, process.lowPtQuadStepTracks)


process.LowPtTripletStepTask = cms.Task(process.lowPtTripletStep, process.lowPtTripletStepClusters, process.lowPtTripletStepHitDoublets, process.lowPtTripletStepHitTriplets, process.lowPtTripletStepSeedLayers, process.lowPtTripletStepSeeds, process.lowPtTripletStepTrackCandidates, process.lowPtTripletStepTrackingRegions, process.lowPtTripletStepTracks)


process.MixedTripletStepTask = cms.Task(process.chargeCut2069Clusters, process.mixedTripletStep, process.mixedTripletStepClassifier1, process.mixedTripletStepClassifier2, process.mixedTripletStepClusters, process.mixedTripletStepHitDoubletsA, process.mixedTripletStepHitDoubletsB, process.mixedTripletStepHitTripletsA, process.mixedTripletStepHitTripletsB, process.mixedTripletStepSeedLayersA, process.mixedTripletStepSeedLayersB, process.mixedTripletStepSeeds, process.mixedTripletStepSeedsA, process.mixedTripletStepSeedsB, process.mixedTripletStepTrackCandidates, process.mixedTripletStepTrackingRegionsA, process.mixedTripletStepTrackingRegionsB, process.mixedTripletStepTracks)


process.PixelLessStepTask = cms.Task(process.pixelLessStep, process.pixelLessStepClassifier1, process.pixelLessStepClassifier2, process.pixelLessStepClusters, process.pixelLessStepHitDoublets, process.pixelLessStepHitTriplets, process.pixelLessStepSeedLayers, process.pixelLessStepSeeds, process.pixelLessStepTrackCandidates, process.pixelLessStepTrackingRegions, process.pixelLessStepTracks)


process.PixelPairStepTask = cms.Task(process.pixelPairStep, process.pixelPairStepClusters, process.pixelPairStepHitDoublets, process.pixelPairStepHitDoubletsB, process.pixelPairStepSeedLayers, process.pixelPairStepSeeds, process.pixelPairStepSeedsA, process.pixelPairStepSeedsB, process.pixelPairStepTrackCandidates, process.pixelPairStepTrackingRegions, process.pixelPairStepTrackingRegionsSeedLayersB, process.pixelPairStepTracks)


process.TobTecStepTask = cms.Task(process.tobTecStep, process.tobTecStepClassifier1, process.tobTecStepClassifier2, process.tobTecStepClusters, process.tobTecStepHitDoubletsPair, process.tobTecStepHitDoubletsTripl, process.tobTecStepHitTripletsTripl, process.tobTecStepSeedLayersPair, process.tobTecStepSeedLayersTripl, process.tobTecStepSeeds, process.tobTecStepSeedsPair, process.tobTecStepSeedsTripl, process.tobTecStepTrackCandidates, process.tobTecStepTrackingRegionsPair, process.tobTecStepTrackingRegionsTripl, process.tobTecStepTracks)


process.beamhaloTracksTask = cms.Task(process.beamhaloTrackCandidates, process.beamhaloTrackerSeedingLayers, process.beamhaloTrackerSeeds, process.beamhaloTracks)


process.caloJetsForTrkTask = cms.Task(process.ak4CaloJetsForTrk, process.caloTowerForTrk)


process.combinatorialcosmicseedinglayersP5Task = cms.Task(process.combinatorialcosmicseedingpairsTECnegP5, process.combinatorialcosmicseedingpairsTECposP5, process.combinatorialcosmicseedingpairsTOBP5, process.combinatorialcosmicseedingtripletsP5)


process.cosmicDCTracksSeqTask = cms.Task(process.cosmicDCCkfTrackCandidates, process.cosmicDCSeeds, process.cosmicDCTracks)


process.cosmictracksP5BottomTask = cms.Task(process.cosmicCandidateFinderP5Bottom, process.cosmicseedfinderP5Bottom, process.cosmictrackfinderP5Bottom)


process.cosmictracksP5Task = cms.Task(process.cosmicCandidateFinderP5, process.cosmicTrackSplitting, process.cosmicseedfinderP5, process.cosmictrackfinderCosmics, process.cosmictrackfinderP5, process.splittedTracksP5)


process.cosmictracksP5TopTask = cms.Task(process.cosmicCandidateFinderP5Top, process.cosmicseedfinderP5Top, process.cosmictrackfinderP5Top)


process.ctfTracksCombinedSeedsTask = cms.Task(process.MixedLayerPairs, process.PixelLayerTriplets, process.ckfTrackCandidatesCombinedSeeds, process.ctfCombinedSeeds, process.globalCombinedSeeds, process.globalSeedsFromPairsWithVertices, process.globalSeedsFromTriplets)


process.ctfTracksNoOverlapsTask = cms.Task(process.ckfTrackCandidatesNoOverlaps, process.ctfNoOverlaps)


process.ctfTracksPixelLessTask = cms.Task(process.ckfTrackCandidatesPixelLess, process.ctfPixelLess, process.globalPixelLessSeeds, process.pixelLessLayerPairs4PixelLessTracking)


process.ctftracksP5BottomTask = cms.Task(process.ckfTrackCandidatesP5Bottom, process.combinatorialcosmicseedfinderP5Bottom, process.combinatorialcosmicseedingpairsTECnegP5Bottom, process.combinatorialcosmicseedingpairsTECposP5Bottom, process.combinatorialcosmicseedingpairsTOBP5Bottom, process.combinatorialcosmicseedingtripletsP5Bottom, process.combinedP5SeedsForCTFBottom, process.ctfWithMaterialTracksP5Bottom, process.simpleCosmicBONSeedingLayersBottom, process.simpleCosmicBONSeedsBottom)


process.ctftracksP5Task = cms.Task(process.ckfTrackCandidatesP5, process.ckfTrackCandidatesP5LHCNavigation, process.combinatorialcosmicseedfinderP5, process.combinatorialcosmicseedinglayersP5Task, process.combinedP5SeedsForCTF, process.ctfWithMaterialTracksCosmics, process.ctfWithMaterialTracksP5, process.ctfWithMaterialTracksP5LHCNavigation, process.simpleCosmicBONSeedingLayers, process.simpleCosmicBONSeeds)


process.ctftracksP5TopTask = cms.Task(process.ckfTrackCandidatesP5Top, process.combinatorialcosmicseedfinderP5Top, process.combinatorialcosmicseedingpairsTECnegP5Top, process.combinatorialcosmicseedingpairsTECposP5Top, process.combinatorialcosmicseedingpairsTOBP5Top, process.combinatorialcosmicseedingtripletsP5Top, process.combinedP5SeedsForCTFTop, process.ctfWithMaterialTracksP5Top, process.simpleCosmicBONSeedingLayersTop, process.simpleCosmicBONSeedsTop)


process.doAlldEdXEstimatorsCTFP5LHCTask = cms.Task(process.dedxHarmonic2CTFP5LHC, process.dedxHitInfoCTFP5LHC, process.dedxTruncated40CTFP5LHC)


process.doAlldEdXEstimatorsCTFTask = cms.Task(process.dedxHarmonic2CTF, process.dedxHitInfoCTF, process.dedxTruncated40CTF)


process.doAlldEdXEstimatorsCosmicTFTask = cms.Task(process.dedxHarmonic2CosmicTF, process.dedxHitInfoCosmicTF, process.dedxTruncated40CosmicTF)


process.doAlldEdXEstimatorsTask = cms.Task(process.dedxHarmonic2, process.dedxHitInfo, process.dedxPixelAndStripHarmonic2T085, process.dedxPixelHarmonic2, process.dedxTruncated40)


process.electronSeedsSeqTask = cms.Task(process.detachedTripletStepSeedClusterMask, process.initialStepSeedClusterMask, process.mixedTripletStepSeedClusterMask, process.newCombinedSeeds, process.pixelLessStepSeedClusterMask, process.pixelPairElectronHitDoublets, process.pixelPairElectronSeedLayers, process.pixelPairElectronSeeds, process.pixelPairElectronTrackingRegions, process.stripPairElectronHitDoublets, process.stripPairElectronSeedLayers, process.stripPairElectronSeeds, process.stripPairElectronTrackingRegions, process.tripletElectronClusterMask, process.tripletElectronHitDoublets, process.tripletElectronHitTriplets, process.tripletElectronSeedLayers, process.tripletElectronSeeds, process.tripletElectronTrackingRegions)


process.generalTracksTask = cms.Task(process.duplicateTrackCandidates, process.duplicateTrackClassifier, process.generalTracks, process.mergedDuplicateTracks)


process.muonSeededStepCoreInOutTask = cms.Task(process.muonSeededSeedsInOut, process.muonSeededTrackCandidatesInOut, process.muonSeededTracksInOut)


process.muonSeededStepCoreOutInTask = cms.Task(process.muonSeededSeedsOutIn, process.muonSeededTrackCandidatesOutIn, process.muonSeededTracksOutIn)


process.muonSeededStepCoreTask = cms.Task(process.muonSeededStepCoreInOutTask, process.muonSeededStepCoreOutInTask)


process.muonSeededStepDebugInOutTask = cms.Task(process.muonSeededSeedsInOutAsTracks, process.muonSeededTrackCandidatesInOutAsTracks)


process.muonSeededStepDebugTask = cms.Task(process.muonSeededSeedsOutInAsTracks, process.muonSeededStepDebugInOutTask, process.muonSeededTrackCandidatesOutInAsTracks)


process.muonSeededStepExtraInOutTask = cms.Task(process.muonSeededTracksInOutClassifier)


process.muonSeededStepExtraTask = cms.Task(process.muonSeededStepExtraInOutTask, process.muonSeededTracksOutInClassifier)


process.muonSeededStepTask = cms.Task(process.earlyMuons, process.muonSeededStepCoreTask, process.muonSeededStepExtraTask)


process.regionalCosmicTracksTask = cms.Task(process.regionalCosmicCkfTrackCandidates, process.regionalCosmicTrackerSeedingLayers, process.regionalCosmicTrackerSeeds, process.regionalCosmicTracks)


process.siPixelClustersPreSplittingTask = cms.Task(process.siPixelClustersPreSplitting)


process.siPixelRecHitsPreSplittingTask = cms.Task(process.siPixelRecHitsPreSplitting)


process.striptrackerlocalrecoTask = cms.Task(process.siStripClusters, process.siStripMatchedRecHits, process.siStripZeroSuppression)


process.trackerlocalrecoBottomTask = cms.Task(process.siPixelClustersBottom, process.siPixelRecHitsBottom, process.siStripClustersBottom, process.siStripMatchedRecHitsBottom, process.topBottomClusterInfoProducerBottom)


process.trackerlocalrecoTopTask = cms.Task(process.siPixelClustersTop, process.siPixelRecHitsTop, process.siStripClustersTop, process.siStripMatchedRecHitsTop, process.topBottomClusterInfoProducerTop)


process.tracksP5BottomTask = cms.Task(process.cosmictracksP5BottomTask, process.ctftracksP5BottomTask)


process.tracksP5TopTask = cms.Task(process.cosmictracksP5TopTask, process.ctftracksP5TopTask)


process.hpsPFTauBasicDiscriminatorsTask = cms.Task(process.hpsPFTauBasicDiscriminators)


process.hpsPFTauBasicDiscriminatorsdR03Task = cms.Task(process.hpsPFTauBasicDiscriminatorsdR03)


process.hpsPFTauDiscriminationByIsolationMVArun2v1DBdR03oldDMwLTTask = cms.Task(process.hpsPFTauDiscriminationByIsolationMVArun2v1DBdR03oldDMwLT, process.hpsPFTauDiscriminationByIsolationMVArun2v1DBdR03oldDMwLTraw)


process.hpsPFTauDiscriminationByIsolationMVArun2v1DBnewDMwLTTask = cms.Task(process.hpsPFTauDiscriminationByIsolationMVArun2v1DBnewDMwLT, process.hpsPFTauDiscriminationByIsolationMVArun2v1DBnewDMwLTraw)


process.hpsPFTauDiscriminationByIsolationMVArun2v1DBoldDMwLTTask = cms.Task(process.hpsPFTauDiscriminationByIsolationMVArun2v1DBoldDMwLT, process.hpsPFTauDiscriminationByIsolationMVArun2v1DBoldDMwLTraw)


process.hpsPFTauMVAIsolation2Task = cms.Task(process.hpsPFTauDiscriminationByIsolationMVArun2v1DBdR03oldDMwLTTask, process.hpsPFTauDiscriminationByIsolationMVArun2v1DBnewDMwLTTask, process.hpsPFTauDiscriminationByIsolationMVArun2v1DBoldDMwLTTask)


process.hpsPFTauVertexAndImpactParametersTask = cms.Task(process.hpsPFTauPrimaryVertexProducer, process.hpsPFTauSecondaryVertexProducer, process.hpsPFTauTransverseImpactParameters)


process.mvaIsolation2TaskRun2 = cms.Task(process.discriminationByIsolationMVArun2v1, process.discriminationByIsolationMVArun2v1raw, process.hpsPFTauBasicDiscriminators)


process.produceHPSPFTausTask = cms.Task(process.hpsPFTauProducer, process.hpsPFTauProducerSansRefs, process.hpsSelectionDiscriminator)


process.recoTauCommonTask = cms.Task(process.ak4PFJetTracksAssociatorAtVertex, process.pfRecoTauTagInfoProducer, process.recoTauAK4PFJets08Region, process.recoTauPileUpVertices)


process.InitialStepTask = cms.Task(process.caloJetsForTrkTask, process.firstStepPrimaryVertices, process.firstStepPrimaryVerticesUnsorted, process.initialStep, process.initialStepClassifier1, process.initialStepHitDoublets, process.initialStepHitQuadruplets, process.initialStepSeedLayers, process.initialStepSeeds, process.initialStepTrackCandidates, process.initialStepTrackCandidatesMkFit, process.initialStepTrackCandidatesMkFitConfig, process.initialStepTrackCandidatesMkFitSeeds, process.initialStepTrackRefsForJets, process.initialStepTrackingRegions, process.initialStepTracks, process.mkFitEventOfHits, process.mkFitGeometryESProducer, process.mkFitSiPixelHits)


process.doAllCosmicdEdXEstimatorsTask = cms.Task(process.doAlldEdXEstimatorsCTFP5LHCTask, process.doAlldEdXEstimatorsCTFTask, process.doAlldEdXEstimatorsCosmicTFTask)


process.iterTrackingEarlyTask = cms.Task(process.DetachedQuadStepTask, process.DetachedTripletStepTask, process.HighPtTripletStepTask, process.InitialStepTask, process.JetCoreRegionalStepTask, process.LowPtQuadStepTask, process.LowPtTripletStepTask, process.MixedTripletStepTask, process.PixelLessStepTask, process.PixelPairStepTask, process.TobTecStepTask)


process.iterTrackingTask = cms.Task(process.ConvStepTask, process.InitialStepPreSplittingTask, process.conversionStepTracks, process.earlyGeneralTracks, process.generalTracksTask, process.iterTrackingEarlyTask, process.muonSeededStepTask, process.preDuplicateMergingGeneralTracks, process.trackdnn_source, process.trackerClusterCheck)


process.pixeltrackerlocalrecoTask = cms.Task(process.siPixelClustersPreSplittingTask, process.siPixelRecHitsPreSplittingTask)


process.trackerCosmics_TopBotTask = cms.Task(process.trackerlocalrecoBottomTask, process.trackerlocalrecoTopTask, process.tracksP5BottomTask, process.tracksP5TopTask)


process.trackerlocalrecoTask = cms.Task(process.clusterSummaryProducer, process.pixeltrackerlocalrecoTask, process.striptrackerlocalrecoTask)


process.tracksP5Task = cms.Task(process.cosmictracksP5Task, process.ctftracksP5Task, process.doAllCosmicdEdXEstimatorsTask, process.siPixelClusterShapeCache)


process.produceAndDiscriminateHPSPFTausTask = cms.Task(process.hpsPFTauBasicDiscriminatorsTask, process.hpsPFTauBasicDiscriminatorsdR03Task, process.hpsPFTauDiscriminationByDeadECALElectronRejection, process.hpsPFTauDiscriminationByDecayModeFinding, process.hpsPFTauDiscriminationByDecayModeFindingNewDMs, process.hpsPFTauDiscriminationByDecayModeFindingOldDMs, process.hpsPFTauDiscriminationByMuonRejection3, process.hpsPFTauVertexAndImpactParametersTask, process.produceHPSPFTausTask)


process.recoTauClassicHPSTask = cms.Task(process.ak4PFJetsLegacyHPSPiZeros, process.ak4PFJetsRecoTauChargedHadrons, process.combinatoricRecoTaus, process.produceAndDiscriminateHPSPFTausTask)


process.ckftracksTask = cms.Task(process.doAlldEdXEstimatorsTask, process.electronSeedsSeqTask, process.iterTrackingTask)


process.ckftracks_plus_pixellessTask = cms.Task(process.ckftracksTask, process.ctfTracksPixelLessTask)


process.ckftracks_woBHTask = cms.Task(process.doAlldEdXEstimatorsTask, process.electronSeedsSeqTask, process.iterTrackingTask)


process.ckftracks_wodEdXTask = cms.Task(process.electronSeedsSeqTask, process.iterTrackingTask)


process.trackingGlobalRecoTask = cms.Task(process.ckftracksTask, process.trackExtrapolator)


process.PFTauTask = cms.Task(process.recoTauClassicHPSTask, process.recoTauCommonTask)


process.Conv2Step = cms.Sequence(process.Conv2StepTask)


process.ConvStep = cms.Sequence(process.ConvStepTask)


process.DetachedQuadStep = cms.Sequence(process.DetachedQuadStepTask)


process.DetachedTripletStep = cms.Sequence(process.DetachedTripletStepTask)


process.DisplacedGeneralStep = cms.Sequence(process.DisplacedGeneralStepTask)


process.HighPtTripletStep = cms.Sequence(process.HighPtTripletStepTask)


process.InitialStep = cms.Sequence(process.InitialStepTask)


process.InitialStepPreSplitting = cms.Sequence(process.InitialStepPreSplittingTask)


process.JetCoreRegionalStep = cms.Sequence(process.JetCoreRegionalStepTask)


process.LowPtQuadStep = cms.Sequence(process.LowPtQuadStepTask)


process.LowPtTripletStep = cms.Sequence(process.LowPtTripletStepTask)


process.MixedTripletStep = cms.Sequence(process.MixedTripletStepTask)


process.PixelLessStep = cms.Sequence(process.PixelLessStepTask)


process.PixelPairStep = cms.Sequence(process.PixelPairStepTask)


process.TobTecStep = cms.Sequence(process.TobTecStepTask)


process.beamhaloTracksSeq = cms.Sequence(process.beamhaloTracksTask)


process.caloJetsForTrk = cms.Sequence(process.caloJetsForTrkTask)


process.ckftracks = cms.Sequence(process.ckftracksTask)


process.ckftracks_plus_pixelless = cms.Sequence(process.ckftracks_plus_pixellessTask)


process.ckftracks_woBH = cms.Sequence(process.ckftracks_woBHTask)


process.ckftracks_wodEdX = cms.Sequence(process.ckftracks_wodEdXTask)


process.combinatorialcosmicseedinglayersP5 = cms.Sequence(process.combinatorialcosmicseedinglayersP5Task)


process.cosmicDCTracksSeq = cms.Sequence(process.cosmicDCTracksSeqTask)


process.cosmictracksP5 = cms.Sequence(process.cosmictracksP5Task)


process.cosmictracksP5Bottom = cms.Sequence(process.cosmictracksP5BottomTask)


process.cosmictracksP5Top = cms.Sequence(process.cosmictracksP5TopTask)


process.ctfTracksCombinedSeeds = cms.Sequence(process.ctfTracksCombinedSeedsTask)


process.ctfTracksNoOverlaps = cms.Sequence(process.ctfTracksNoOverlapsTask)


process.ctfTracksPixelLess = cms.Sequence(process.ctfTracksPixelLessTask)


process.ctftracksP5 = cms.Sequence(process.ctftracksP5Task)


process.ctftracksP5Bottom = cms.Sequence(process.ctftracksP5BottomTask)


process.ctftracksP5Top = cms.Sequence(process.ctftracksP5TopTask)


process.doAllCosmicdEdXEstimators = cms.Sequence(process.doAllCosmicdEdXEstimatorsTask)


process.doAlldEdXEstimators = cms.Sequence(process.doAlldEdXEstimatorsTask)


process.doAlldEdXEstimatorsCTF = cms.Sequence(process.doAlldEdXEstimatorsCTFTask)


process.doAlldEdXEstimatorsCTFP5LHC = cms.Sequence(process.doAlldEdXEstimatorsCTFP5LHCTask)


process.doAlldEdXEstimatorsCosmicTF = cms.Sequence(process.doAlldEdXEstimatorsCosmicTFTask)


process.electronSeedsSeq = cms.Sequence(process.electronSeedsSeqTask)


process.generalTracksSequence = cms.Sequence(process.generalTracksTask)


process.iterTracking = cms.Sequence(process.iterTrackingTask)


process.iterTrackingEarly = cms.Sequence(process.iterTrackingEarlyTask)


process.muonSeededStep = cms.Sequence(process.muonSeededStepTask)


process.muonSeededStepCore = cms.Sequence(process.muonSeededStepCoreTask)


process.muonSeededStepCoreInOut = cms.Sequence(process.muonSeededStepCoreInOutTask)


process.muonSeededStepCoreOutIn = cms.Sequence(process.muonSeededStepCoreOutInTask)


process.muonSeededStepDebug = cms.Sequence(process.muonSeededStepDebugTask)


process.muonSeededStepDebugInOut = cms.Sequence(process.muonSeededStepDebugInOutTask)


process.muonSeededStepExtra = cms.Sequence(process.muonSeededStepExtraTask)


process.muonSeededStepExtraInOut = cms.Sequence(process.muonSeededStepExtraInOutTask)


process.pixeltrackerlocalreco = cms.Sequence(process.pixeltrackerlocalrecoTask)


process.regionalCosmicTracksSeq = cms.Sequence(process.regionalCosmicTracksTask)


process.striptrackerlocalreco = cms.Sequence(process.striptrackerlocalrecoTask)


process.trackerCosmics_TopBot = cms.Sequence(process.trackerCosmics_TopBotTask)


process.trackerlocalreco = cms.Sequence(process.trackerlocalrecoTask)


process.trackerlocalrecoBottom = cms.Sequence(process.trackerlocalrecoBottomTask)


process.trackerlocalrecoTop = cms.Sequence(process.trackerlocalrecoTopTask)


process.trackingGlobalReco = cms.Sequence(process.trackingGlobalRecoTask)


process.tracksP5 = cms.Sequence(process.tracksP5Task)


process.tracksP5Bottom = cms.Sequence(process.tracksP5BottomTask)


process.tracksP5Top = cms.Sequence(process.tracksP5TopTask)


process.tracksP5_wodEdX = cms.Sequence(process.tracksP5Task)


process.PFTau = cms.Sequence(process.PFTauTask)


process.hpsPFTauMVAIsolation2Seq = cms.Sequence(process.hpsPFTauBasicDiscriminatorsTask, process.hpsPFTauBasicDiscriminatorsdR03Task, process.hpsPFTauMVAIsolation2Task)


process.hpsPFTauVertexAndImpactParametersSeq = cms.Sequence(process.hpsPFTauVertexAndImpactParametersTask)


process.mvaIsolation2SeqRun2 = cms.Sequence(process.mvaIsolation2TaskRun2)


process.produceAndDiscriminateHPSPFTaus = cms.Sequence(process.produceAndDiscriminateHPSPFTausTask)


process.produceHPSPFTaus = cms.Sequence(process.produceHPSPFTausTask)


process.recoTauClassicHPSSequence = cms.Sequence(process.recoTauClassicHPSTask)


process.recoTauCommonSequence = cms.Sequence(process.recoTauCommonTask)


process.p = cms.Path(process.siStripMatchedRecHits+process.siPixelRecHits+process.MeasurementTrackerEvent+process.TrackRefitter+process.hltFilter+process.fevt)


