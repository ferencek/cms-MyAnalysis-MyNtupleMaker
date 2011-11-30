import FWCore.ParameterSet.Config as cms

GenParticles = cms.EDProducer("MyNtupleMaker_GenParticles",
    InputTag = cms.InputTag('genParticles'),
    Prefix = cms.string(''),
    Suffix = cms.string(''),
    MaxSize = cms.int32(20), # turned off if negative
    PdgIDsOfInterest = cms.vint32( # will be stored regardless of the MaxSize parameter
        # Full list of PDG ID codes can be found in http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/SimGeneral/HepPDTESSource/data/pythiaparticle.tbl
        4,    # c
        5,    # b
        13,   # mu
        411,  # D+
        421,  # D0
        431,  # D_s+
        4122, # Lambda_c+
        4132, # Xi_c0
        4232, # Xi_c+
        4332, # Omega_c0
        511,  # B0
        521,  # B+
        531,  # B_s0
        541,  # B_c+
        5122, # Lambda_b0
        5132, # Xi_b-
        5232, # Xi_b0
        5332, # Omega_b-
    )
)
