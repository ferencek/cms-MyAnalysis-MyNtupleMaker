###############################
####### Parameters ############
###############################
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('python')

options.register('runOnData',
    False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)

options.register('globalTag',
    'START42_V14B::All',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Global tag to be used"
)

options.register('outputFilename',
    'MyNtupleMaker_output.root',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Name of the output file"
)

options.register('outputPrefix',
    '',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Prefix for the output file names"
)

options.register('produceSkim',
    False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Switch to turn ON/OFF skim production"
)

options.register('reportEvery',
    10,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Report every N events (default is N=10)"
)
## 'maxEvents' is already registered by the Framework, changing default value
options.setDefault('maxEvents', 10)

options.parseArguments()

print "Run on data: %s"%('True' if options.runOnData else 'False')

## For debugging
#print options

## Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *

from PhysicsTools.PatAlgos.tools.coreTools import *

## Remove certain objects from the default sequence
removeAllPATObjectsBut(process, ['Jets', 'Muons', 'METs'])

############## IMPORTANT ########################################
# If you run over many samples and you save the log, remember to reduce
# the size of the output by prescaling the report of the event number
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
process.MessageLogger.cerr.default.limit = 10
#################################################################

## JEC levels
inputJetCorrLevelsCalo = cms.vstring(['L1Offset', 'L2Relative', 'L3Absolute'])
inputJetCorrLevelsPF   = cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute'])

## Load MyNtupleMaker modules
process.load('MyAnalysis.MyNtupleMaker.MyNtupleMaker_cff')

## Make sure a correct global tag is used (please refer to https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions#Valid_Global_Tags_by_Release)
process.GlobalTag.globaltag = options.globalTag

## Events to process
process.maxEvents.input = options.maxEvents

## Options and Output Report
process.options.wantSummary = cms.untracked.bool(True)
#process.options.makeTriggerResults = cms.untracked.bool(True)

## Input files
process.source.fileNames = [
    '/store/mc/Summer11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/02719D6B-1398-E011-AA71-001A92971B94.root'
]

## If running on data
if options.runOnData:
    ## Turn off MC matching for the process
    removeMCMatching(process, ['All'])

    inputJetCorrLevelsCalo = cms.vstring(['L1Offset', 'L2Relative', 'L3Absolute', 'L2L3Residual'])
    inputJetCorrLevelsPF   = cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'])
 
    process.source.fileNames = [
        '/store/data/Run2011A/Jet/AOD/May10ReReco-v1/0000/94A6E942-447C-E011-8F5F-0024E8768D68.root'
    ]
    ## Get new calibration for the jet probability b-tagger (already included in START42_V14B::All hence only needed for data)
    ## The following modifications are based on https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagJetProbabilityCalibration#Use_the_new_calibration_in_42x_r
    process.GlobalTag.toGet = cms.VPSet(
        cms.PSet(
            record = cms.string("BTagTrackProbability2DRcd"),
            tag = cms.string("TrackProbabilityCalibration_2D_2011Data_v1_offline"),
            connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_BTAU")
        ),
        cms.PSet(
            record = cms.string("BTagTrackProbability3DRcd"),
            tag = cms.string("TrackProbabilityCalibration_3D_2011Data_v1_offline"),
            connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_BTAU")
        )
    )

## For L1FastJet Corrections ( https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections )
##-------------------- Import the Jet RECO modules -----------------------
process.load('RecoJets.Configuration.RecoPFJets_cff')
##-------------------- Turn-on the FastJet density calculation -----------------------
#process.kt6PFJets.doRhoFastjet = True # added automatically by PAT when L1FastJet is used
##-------------------- Turn-on the FastJet jet area calculation for your favorite algorithm -----------------------
process.ak5PFJets.doAreaFastjet = True
process.ak7PFJets.doAreaFastjet = True

## Calculate the FastJet desity to correct the isolation
process.kt6PFJetsForIsolation = process.kt6PFJets.clone()
process.kt6PFJetsForIsolation.doRhoFastjet = True
process.kt6PFJetsForIsolation.Rho_EtaMax = 2.5

from PhysicsTools.PatAlgos.tools.jetTools import *
## The default jet collection is ak5 CaloJets
switchJetCollection(process,cms.InputTag('ak5CaloJets'),
    doJTA        = True,
    doBTagging   = True,
    jetCorrLabel = ('AK5Calo', inputJetCorrLevelsCalo),
    doType1MET   = False,
    genJetCollection=cms.InputTag("ak5GenJets"),
    doJetID      = True
)
## Add ak7 CaloJets
addJetCollection(process,cms.InputTag('ak7CaloJets'),
    'AK7', 'Calo',
    doJTA        = True,
    doBTagging   = True,
    jetCorrLabel = ('AK7Calo', inputJetCorrLevelsCalo),
    doType1MET   = False,
    doL1Cleaning = False,
    doL1Counters = False,
    genJetCollection = cms.InputTag("ak7GenJets"),
    doJetID      = True,
    jetIdLabel   = 'ak7'
)
## Add ak5 PFJets
addJetCollection(process,cms.InputTag('ak5PFJets'),
    'AK5', 'PF',
    doJTA        = True,
    doBTagging   = True,
    jetCorrLabel = ('AK5PF', inputJetCorrLevelsPF),
    doType1MET   = False,
    doL1Cleaning = False,
    doL1Counters = False,
    genJetCollection=cms.InputTag("ak5GenJets"),
    doJetID      = False
)
## Add ak7 PFJets
addJetCollection(process,cms.InputTag('ak7PFJets'),
    'AK7', 'PF',
    doJTA        = True,
    doBTagging   = True,
    jetCorrLabel = ('AK7PF', inputJetCorrLevelsPF),
    doType1MET   = False,
    doL1Cleaning = False,
    doL1Counters = False,
    genJetCollection = cms.InputTag("ak7GenJets"),
    doJetID      = False
)

## Define jet selection
process.selectedPatJets.cut        = 'pt > 10.0'
process.selectedPatJetsAK5PF.cut   = 'pt > 10.0'
process.selectedPatJetsAK7Calo.cut = 'pt > 10.0'
process.selectedPatJetsAK7PF.cut   = 'pt > 10.0'

## Add PFMET
from PhysicsTools.PatAlgos.tools.metTools import *
addPfMET(process, 'PF')

## Read JEC uncertainties (might not be available in some global tags)
process.AK5PFJets.ReadJECUncertainty   = True
process.AK5CaloJets.ReadJECUncertainty = True
process.AK7PFJets.ReadJECUncertainty   = True
process.AK7CaloJets.ReadJECUncertainty = True

## Load HBHENoiseFilterResultProducer
process.load('CommonTools/RecoAlgos/HBHENoiseFilterResultProducer_cfi')
# Check the latest recommendation from https://twiki.cern.ch/twiki/bin/view/CMS/HBHEAnomalousSignals2011
process.HBHENoiseFilterResultProducer.minRatio = cms.double(-999)
process.HBHENoiseFilterResultProducer.maxRatio = cms.double(999)
process.HBHENoiseFilterResultProducer.minHPDHits = cms.int32(17)
process.HBHENoiseFilterResultProducer.minRBXHits = cms.int32(999)
process.HBHENoiseFilterResultProducer.minHPDNoOtherHits = cms.int32(10)
process.HBHENoiseFilterResultProducer.minZeros = cms.int32(10)
process.HBHENoiseFilterResultProducer.minHighEHitTime = cms.double(-9999.0)
process.HBHENoiseFilterResultProducer.maxHighEHitTime = cms.double(9999.0)
process.HBHENoiseFilterResultProducer.maxRBXEMF = cms.double(-999.0)
process.HBHENoiseFilterResultProducer.minNumIsolatedNoiseChannels = cms.int32(9999)
process.HBHENoiseFilterResultProducer.minIsolatedNoiseSumE = cms.double(9999)
process.HBHENoiseFilterResultProducer.minIsolatedNoiseSumEt = cms.double(9999)
process.HBHENoiseFilterResultProducer.useTS4TS5 = cms.bool(True)

## EcalDeadCell simpleDeltaR filters
process.load('JetMETAnalysis.simpleDRfilter.simpleDRfilter_cfi')
process.simpleDRfilter.debug = cms.untracked.bool(False)
process.simpleDRfilter.jetInputTag = cms.InputTag("ak5PFJets")
process.simpleDRfilter.metInputTag = cms.InputTag("pfMet")
process.simpleDRfilter.doFilter = cms.untracked.bool(False) # to enable filter
process.simpleDRfilter.makeProfileRoot = False

## Good vertex selection
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
  vertexCollection = cms.InputTag('offlinePrimaryVertices'),
  minimumNDOF      = cms.uint32  (4), # this is > 4
  maxAbsZ          = cms.double  (24.0),
  maxd0            = cms.double  (2.0),
)

## Removal of beam scraping events
process.scrapingVeto = cms.EDFilter("FilterOutScraping",
  applyfilter = cms.untracked.bool  (True),
  debugOn     = cms.untracked.bool  (False),
  numtrack    = cms.untracked.uint32(10),
  thresh      = cms.untracked.double(0.25)
)

## Produce a collection of good primary vertices
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
process.goodOfflinePrimaryVertices = cms.EDFilter("PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone(
        minNdof = cms.double(4.0), # this is >= 4
        maxZ = cms.double(24.0),
        maxRho = cms.double(2.0)
    ),
    src = cms.InputTag('offlinePrimaryVertices')
)

## Event counters
process.nEventsTotal = cms.EDProducer("EventCountProducer")

## MyAnalyzer configuration
process.myAnalyzer = cms.EDFilter('MyAnalyzer',
    fillAllSameLevelAndLowerLevelCuts = cms.untracked.bool(False), # to disable automatic creation of less frequently used histograms
    fillAllCuts                       = cms.untracked.bool(False), # to disable automatic creation of less frequently used histograms
    skimMode                          = cms.untracked.bool(options.produceSkim), # when enabled, the cutEfficiency file and output histograms are not produced
    HLTInputTag                       = cms.InputTag('TriggerResults','','HLT'),
    skimWasMade                       = cms.bool(True),
    eventCounterInputTag              = cms.untracked.InputTag('nEventsTotal'),
    inputCutFile                      = cms.string('cutFile.txt'),
    outputCutEfficiencyFile           = cms.string(((options.outputPrefix + '__') if options.outputPrefix != '' else '') + 'cutEfficiency.txt')
)

## Path definition
process.p = cms.Path(
    process.nEventsTotal*
    process.primaryVertexFilter*
    process.scrapingVeto*
    (
    process.goodOfflinePrimaryVertices+
    process.HBHENoiseFilterResultProducer+
    process.kt6PFJetsForIsolation+
    process.ak5PFJets+
    process.ak7PFJets+
    process.simpleDRfilter
    )*
    process.patDefaultSequence*
    (
    process.AK5CaloJets+
    process.AK7CaloJets+
    process.AK5GenJets+
    process.AK7GenJets+
    process.AK5PFJets+
    process.AK7PFJets+
    process.CaloMET+
    process.EventSelection+
    process.GenEventInfo+
    process.GenParticles+
    process.Muons+
    process.PFMET+
    process.Vertices
    )*
    process.myAnalyzer
)

## Output file
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string(options.outputFilename),
    # save only events passing the full path
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    dropMetaData = cms.untracked.string("ALL"),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_gtDigis_*_*',
        'keep *_TriggerResults_*_*',
        'drop *_TriggerResults_*_RECO',
        'drop *_TriggerResults_*_PAT',
        'keep *_hltTriggerSummaryAOD_*_*',
        'keep *_nEventsTotal_*_*',
        'keep *_kt6PFJetsAK5PF_rho_*',
        'keep *_kt6PFJetsForIsolation_rho_*',
        'keep *_AK5CaloJets_*_*',
        'keep *_AK7CaloJets_*_*',
        'keep *_AK5GenJets_*_*',
        'keep *_AK7GenJets_*_*',
        'keep *_AK5PFJets_*_*',
        'keep *_AK7PFJets_*_*',
        'keep *_CaloMET_*_*',
        'keep *_EventSelection_*_*',
        'keep *_GenEventInfo_*_*',
        'keep *_GenParticles_*_*',
        'keep *_Muons_*_*',
        'keep *_PFMET_*_*',
        'keep *_Vertices_*_*'
    )
)

## EndPath definition
process.outpath = cms.EndPath(process.out)

## Schedule definition
process.schedule = cms.Schedule(process.p)
process.schedule.append(process.outpath)
