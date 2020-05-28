import FWCore.ParameterSet.Config as cms



process = cms.Process("Demo")



### standard MessageLoggerConfiguration
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use#
    fileNames = cms.untracked.vstring(
'/store/mc/RunIISummer16MiniAODv3/JpsiToMuMu_JpsiPt8_TuneCUEP8M1_13TeV-pythia8/MINIAODSIM/NoPU_94X_mcRun2_asymptotic_v3_ext1-v2/60000/52CC120C-1300-EA11-A8EE-008CFAC91A30.root'
    )
)

process.analysis = cms.EDAnalyzer('MuonCalibAnalyzer',
                                  muons = cms.InputTag("slimmedMuons"),
                                  met = cms.InputTag("slimmedMETs"),
                                  genParticles = cms.InputTag("prunedGenParticles"),
                                  vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                  isOnia = cms.bool(True)
)


process.p = cms.Path(process.analysis)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')


#import PhysicsTools.PythonAnalysis.LumiList as LumiList
#jsonFile='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Reprocessing/Cert_13TeV_27Jan2016ReReco_Collisions15_ZeroTesla_25ns_JSON_MuonPhys.txt'
#myLumis = LumiList.LumiList(filename = jsonFile).getCMSSWString().split(',')
#process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange()
#process.source.lumisToProcess.extend(myLumis)


