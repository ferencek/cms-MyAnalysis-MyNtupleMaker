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
    'START42_V13::All',
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
removeAllPATObjectsBut(process, ['Jets', 'Muons'])

############## IMPORTANT ########################################
# If you run over many samples and you save the log, remember to reduce
# the size of the output by prescaling the report of the event number
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
process.MessageLogger.cerr.default.limit = 10
#################################################################

## JEC levels
inputJetCorrLevelsCalo = cms.vstring(['L1Offset', 'L2Relative', 'L3Absolute'])
inputJetCorrLevelsPF   = cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute'])

if options.runOnData:
    inputJetCorrLevelsCalo = cms.vstring(['L1Offset', 'L2Relative', 'L3Absolute', 'L2L3Residual'])
    inputJetCorrLevelsPF   = cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'])

## Load MyNtupleMaker modules
process.load('MyAnalysis.MyNtupleMaker.MyNtupleMaker_cff')

## Make sure a correct global tag is used (please refer to https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions#Valid_Global_Tags_by_Release)
process.GlobalTag.globaltag = options.globalTag
#process.GlobalTag.globaltag = 'START42_V13::All' # ===> First complete JEC set for 42x 2011 data (https://indico.cern.ch/getFile.py/access?contribId=8&resId=0&materialId=slides&confId=143981)
#process.GlobalTag.globaltag = 'START42_V12::All' # ===> for Summer11 MC analyzed in 42X (contains Jec11_V1, does not contain "residual" JEC and uncertainties yet...)
#process.GlobalTag.globaltag = 'START41_V0::All' # ===> for 41X MC analyzed in 41X (contains Jec10V3)

## Events to process
process.maxEvents.input = options.maxEvents

## Options and Output Report
process.options.wantSummary = cms.untracked.bool(True)
#process.options.makeTriggerResults = cms.untracked.bool(True)

## Input files
process.source.fileNames = [
    '/store/mc/Summer11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/02719D6B-1398-E011-AA71-001A92971B94.root'
]

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

## Read JEC uncertainties (might not be available in some global tag)
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

# Good vertex selection
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
  vertexCollection = cms.InputTag('offlinePrimaryVertices'),
  minimumNDOF      = cms.uint32  (4),
  maxAbsZ          = cms.double  (24.0),
  maxd0            = cms.double  (2.0),
)

# Removal of beam scraping events
process.scrapingVeto = cms.EDFilter("FilterOutScraping",
  applyfilter = cms.untracked.bool  (False),
  debugOn     = cms.untracked.bool  (False),
  numtrack    = cms.untracked.uint32(10),
  thresh      = cms.untracked.double(0.25)
)

## Event counters
process.nEventsTotal = cms.EDProducer("EventCountProducer")

## Path definition
process.p = cms.Path(
    process.nEventsTotal*
    process.primaryVertexFilter*
    process.scrapingVeto*
    process.HBHENoiseFilterResultProducer*
    #process.kt6PFJets* # added automatically by PAT when L1FastJet is used
    process.kt6PFJetsForIsolation*
    process.ak5PFJets*
    process.ak7PFJets*
    process.simpleDRfilter*
    process.patDefaultSequence*
    (
    process.AK5CaloJets+
    process.AK7CaloJets+
    process.AK5GenJets+
    process.AK7GenJets+
    process.AK5PFJets+
    process.AK7PFJets+
    process.EventSelection+
    process.GenEventInfo+
    process.GenParticles+
    process.Muons+
    process.Vertices
    )
)

## Output file
process.output = cms.OutputModule("PoolOutputModule",
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
        'keep *_nEventsTotal_*_*',
        'keep *_kt6PFJets_rho_*',
        'keep *_kt6PFJetsForIsolation_rho_*',
        'keep *_AK5CaloJets_*_*',
        'keep *_AK7CaloJets_*_*',
        'keep *_AK5GenJets_*_*',
        'keep *_AK7GenJets_*_*',
        'keep *_AK5PFJets_*_*',
        'keep *_AK7PFJets_*_*',
        'keep *_EventSelection_*_*',
        'keep *_GenEventInfo_*_*',
        'keep *_GenParticles_*_*',
        'keep *_Muons_*_*',
        'keep *_Vertices_*_*'
    )
)

## Delete predefined output module and Endpath (needed for running with CRAB)
del process.out
del process.outpath

## EndPath definition
process.outpath = cms.EndPath(process.output)

## Schedule definition
process.schedule = cms.Schedule(process.p)
process.schedule.append(process.outpath)
