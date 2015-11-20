####################################################################################################
#MCInfo
####################################################################################################
#cross-section: the cross-section of the MC process in pb-1, 
#twiki: https://twiki.cern.ch/twiki/bin/view/CMS/StandardModelCrossSectionsat13TeV
####################################################################################################
#Number of processed events: Number of events processed by CRAB. Not identical to number of produced events if not all grid jobs succeed!
####################################################################################################
#Number of produced events: Number of total events in the parent sample, i.e. 
#https://cmsweb.cern.ch/das/request?view=list&limit=10&instance=cms_dbs_prod_global&input=dataset%3D%2FElectronHad%2FRun2011A-08Nov2011-v1%2FAOD
####################################################################################################
#Number of selected events: Number of events passing the pre-selection at the nTuple process
####################################################################################################
#Naming of the samples:
#Please use the same naming as in (names array)
#https://svnweb.cern.ch/trac/bat/browser/trunk/AnalysisTools/interface/DataTypes.h
####################################################################################################

datasetInfo = {}

#Data samples
datasetInfo['ElectronHad'] = {"cross-section": 0, "NumberOfProcessedEvents":0}
# datasetInfo['MuHad'] = {"cross-section": 0, "NumberOfProcessedEvents":0}

datasetInfo['SingleElectron'] = {"cross-section": 0, "NumberOfProcessedEvents":0}
datasetInfo['SingleMuon'] = {"cross-section": 0, "NumberOfProcessedEvents":0}

# datasetInfo['SingleElectron_PromptReco_RunD'] = {"cross-section": 0, "NumberOfProcessedEvents":0}
# datasetInfo['SingleElectron_ReReco_RunD'] = {"cross-section": 0, "NumberOfProcessedEvents":0}
# datasetInfo['SingleElectron_ReReco_RunC'] = {"cross-section": 0, "NumberOfProcessedEvents":0}
# datasetInfo['SingleMuon_PromptReco_RunD'] = {"cross-section": 0, "NumberOfProcessedEvents":0}
# datasetInfo['SingleMuon_ReReco_RunD'] = {"cross-section": 0, "NumberOfProcessedEvents":0}
# datasetInfo['SingleMuon_ReReco_RunC'] = {"cross-section": 0, "NumberOfProcessedEvents":0}

# # Cross sections from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO

datasetInfo['TTJets_PowhegPythia8'] = { "cross-section" : 831.76, "NumberOfProcessedEvents" : 19757190 }
datasetInfo['TTJets_PowhegPythia8_scaledown'] = { "cross-section" : 831.76, "NumberOfProcessedEvents" : 9860774 }
datasetInfo['TTJets_PowhegPythia8_scaleup'] = { "cross-section" : 831.76, "NumberOfProcessedEvents" : 9921174 }
datasetInfo['TTJets_PowhegPythia8_mtop1695'] = { "cross-section" : 831.76, "NumberOfProcessedEvents" : 9969800 }
datasetInfo['TTJets_PowhegPythia8_mtop1755'] = { "cross-section" : 831.76, "NumberOfProcessedEvents" : 9904200 }


datasetInfo['TTJets_amcatnloFXFX'] = { "cross-section" : 831.76, "NumberOfProcessedEvents" : 14188545 }
datasetInfo['TTJets_amcatnloFXFX_scaledown'] = { "cross-section" : 831.76, "NumberOfProcessedEvents" : 4986572 }#not available
datasetInfo['TTJets_amcatnloFXFX_scaleup'] = { "cross-section" : 831.76, "NumberOfProcessedEvents" : 4998339 }#not available
datasetInfo['TTJets_amcatnloFXFX_mtop1695'] = { "cross-section" : 831.76, "NumberOfProcessedEvents" : 3203436 }
datasetInfo['TTJets_amcatnloFXFX_mtop1755'] = { "cross-section" : 831.76, "NumberOfProcessedEvents" : 3638767 }

datasetInfo['TTJets_madgraphMLM'] = { "cross-section" : 831.76, "NumberOfProcessedEvents" : 11344206 }
# datasetInfo['TTJets_powhegPythia6'] = { "cross-section" : 831.76, "NumberOfProcessedEvents" : 1913779 }
# datasetInfo['TTJets_PowhegHerwigpp'] = { "cross-section" : 831.76, "NumberOfProcessedEvents" : 4999975 }
datasetInfo['TTJets_amcatnloHerwigpp'] = { "cross-section" : 831.76, "NumberOfProcessedEvents" : 11110130 }

# # Cross sections from https://twiki.cern.ch/twiki/bin/view/CMS/StandardModelCrossSectionsat13TeV
datasetInfo['DYJetsToLL_M50'] = { "cross-section" : 6025.2, "NumberOfProcessedEvents" : 18969962 }
datasetInfo['WJetsToLNu'] = { "cross-section" : 61526.7, "NumberOfProcessedEvents" : 16541248 }

# # Cross sections from https://twiki.cern.ch/twiki/bin/viewauth/CMS/SingleTopSigma
# datasetInfo['TToLeptons_s'] = { "cross-section" : 1, "NumberOfProcessedEvents" : 1 }
# datasetInfo['TBarToLeptons_s'] = { "cross-section" : 1, "NumberOfProcessedEvents" : 1 }
# datasetInfo['TBarToLeptons_t'] = { "cross-section" : 1, "NumberOfProcessedEvents" : 1 }

datasetInfo['TToLeptons_t'] = { "cross-section" : 136.02, "NumberOfProcessedEvents" : 1119438 }
datasetInfo['T_tw'] = { "cross-section" : 35.6, "NumberOfProcessedEvents" : 995600 }
datasetInfo['Tbar_tw'] = { "cross-section" : 35.6, "NumberOfProcessedEvents" : 988500 }

# # G+Jet
datasetInfo['GJets_15toInf'] = { "cross-section" : 365896, "NumberOfProcessedEvents" : 1972730 }

# # QCD EM-enriched
# # Cross sections from https://cms-pdmv.cern.ch/mcm/requests?dataset_name=QCD_*EMEnriched*TuneCUETP8M1*13*&page=0&shown=262173&limit=100
datasetInfo['QCD_EMEnriched_15to20'] = { "cross-section" : 1273000000 * 0.0002, "NumberOfProcessedEvents" : 2122959 }
datasetInfo['QCD_EMEnriched_20to30'] = { "cross-section" : 557600000 * 0.0096, "NumberOfProcessedEvents" : 9184197 }
datasetInfo['QCD_EMEnriched_30to50'] = { "cross-section" : 136000000 * 0.073, "NumberOfProcessedEvents" : 4736965 }
datasetInfo['QCD_EMEnriched_50to80'] = { "cross-section" : 19800000 * 0.146, "NumberOfProcessedEvents" : 5291543 }
datasetInfo['QCD_EMEnriched_80to120'] = { "cross-section" : 2800000 * 0.125, "NumberOfProcessedEvents" : 8130424 }
datasetInfo['QCD_EMEnriched_120to170'] = { "cross-section" : 477000 * 0.132, "NumberOfProcessedEvents" : 8513976 }
datasetInfo['QCD_EMEnriched_170to300'] = { "cross-section" : 114000 * 0.165, "NumberOfProcessedEvents" : 5735584 }
datasetInfo['QCD_EMEnriched_300toInf'] = { "cross-section" : 9000 * 0.15, "NumberOfProcessedEvents" : 3707833 }
# # QCD BCToE
# # Cross sections from https://cms-pdmv.cern.ch/mcm/requests?dataset_name=QCD_Pt_*_bcToE_TuneCUETP8M1_*13*&page=0&shown=262173&limit=100
datasetInfo['QCD_bcToE_15to20'] = { "cross-section" : 1272980000 * 0.0002, "NumberOfProcessedEvents" : 907678 }
datasetInfo['QCD_bcToE_20to30'] = { "cross-section" : 557627000 * 0.00059, "NumberOfProcessedEvents" : 3252331 }
datasetInfo['QCD_bcToE_30to80'] = { "cross-section" : 159068000 * 0.00255, "NumberOfProcessedEvents" : 5097805 }
datasetInfo['QCD_bcToE_80to170'] = { "cross-section" : 3221000 * 0.01183, "NumberOfProcessedEvents" : 5141746 }
datasetInfo['QCD_bcToE_170to250'] = { "cross-section" : 105771 * 0.02492, "NumberOfProcessedEvents" : 2882785 }
datasetInfo['QCD_bcToE_250toInf'] = { "cross-section" : 21094.1 * 0.03375, "NumberOfProcessedEvents" : 2951911 }
# # QCD Mu-enriched
# # Cross sections from https://cms-pdmv.cern.ch/mcm/requests?dataset_name=QCD_Pt*_MuEnrichedPt5_TuneCUETP8M1_*13*&page=0&shown=262173&limit=100
datasetInfo['QCD_MuEnriched_15to20'] = { "cross-section" : 1273190000 * 0.003, "NumberOfProcessedEvents" : 2362016 }
datasetInfo['QCD_MuEnriched_20to30'] = { "cross-section" : 558528000 * 0.0053, "NumberOfProcessedEvents" : 5311955 }
datasetInfo['QCD_MuEnriched_30to50'] = { "cross-section" : 139803000 * 0.01182, "NumberOfProcessedEvents" : 4946084 }
datasetInfo['QCD_MuEnriched_50to80'] = { "cross-section" : 19222500 * 0.02276, "NumberOfProcessedEvents" : 5058536 }
datasetInfo['QCD_MuEnriched_80to120'] = { "cross-section" : 2758420 * 0.03844, "NumberOfProcessedEvents" : 3882463 }
datasetInfo['QCD_MuEnriched_120to170'] = { "cross-section" : 469797 * 0.05362, "NumberOfProcessedEvents" : 4026104 }
datasetInfo['QCD_MuEnriched_170to300'] = { "cross-section" : 117989 * 0.07335, "NumberOfProcessedEvents" : 3949716 }
datasetInfo['QCD_MuEnriched_300to470'] = { "cross-section" : 7820.25 * 0.10196, "NumberOfProcessedEvents" : 3960205 }
datasetInfo['QCD_MuEnriched_470to600'] = { "cross-section" : 645.528 * 0.12242, "NumberOfProcessedEvents" : 1928595 }
datasetInfo['QCD_MuEnriched_600to800'] = { "cross-section" : 187.109 * 0.13412, "NumberOfProcessedEvents" : 1 } #not available 1982667
datasetInfo['QCD_MuEnriched_800to1000'] = { "cross-section" : 32.3486 * 0.14552, "NumberOfProcessedEvents" : 1982314 }
datasetInfo['QCD_MuEnriched_1000toInf'] = { "cross-section" : 10.4305 * 0.15544, "NumberOfProcessedEvents" : 1981954 }
