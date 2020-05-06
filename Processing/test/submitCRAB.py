import os


samples={}

samples['jpsi']=[

'/JPsiToMuMu_Pt20to100-pythia8-gun/RunIISummer16MiniAODv2-PUMoriond17RAW_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',
'/JpsiToMuMu_JpsiPt8_TuneCUEP8M1_13TeV-pythia8-photos/RunIISummer16MiniAODv3-NoPU_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM',
'/JpsiToMuMu_JpsiPt8_TuneCUEP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',
'/JpsiToMuMu_JpsiPt8_TuneCUEP8M1_13TeV-pythia8/RunIISummer16MiniAODv3-NoPU_94X_mcRun2_asymptotic_v3_ext1-v1/MINIAODSIM',
'/JpsiToMuMu_JpsiPt8_TuneCUEP8M1_13TeV-pythia8/RunIISummer16MiniAODv3-NoPU_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM'

]

samples['z']=[
'/ZJToMuMu_mWPilot_TuneCP5_13TeV-powheg-MiNNLO-pythia8-photos/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM',
'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM'
]
json='https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt'

for resonance,slist in samples.iteritems():
    for i,s in enumerate(slist):
        filename="crabSubmit_"+resonance+"_"+str(i)+'.py'
        f=open(filename,'w')
        cfg="""
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()


config.General.requestName = '{name}'
config.General.workArea = 'CRAB'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '{cfgFile}'
config.JobType.outputFiles = ['muonTree.root']
config.Data.inputDataset = '{dataset}'
config.Data.inputDBS = 'global'

config.Data.splitting = 'FileBased'
config.Data.allowNonValidInputDataset = True
config.Data.unitsPerJob = 5
#config.Data.lumiMask = '{json}'
config.Data.outLFNDirBase = '/store/user/emanca/KaMuCa/'
config.Data.publication = False
config.Site.storageSite = 'T2_IT_Pisa'

""".format( name=resonance+"_"+str(i),cfgFile='runZ.py' if resonance=='z' else 'runOnia.py',dataset=s,json=json)

        f.write(cfg)
        f.close()



        os.system('crab submit --config='+filename)
