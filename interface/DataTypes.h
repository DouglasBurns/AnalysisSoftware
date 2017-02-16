/*
 * Enumerators.h
 *
 *  Created on: Jun 30, 2010
 *      Author: lkreczko
 */

#ifndef ENUMERATORS_H_
#define ENUMERATORS_H_
#include <string>
#include <boost/array.hpp>

namespace BAT {
namespace DataType {
enum value {
	ElectronHad,
	// MuHad,
	SingleElectron_RunB,
	SingleElectron_RunC,
	SingleElectron_RunD,
	SingleElectron_RunE,
	SingleElectron_RunF,
	SingleElectron_RunG,
	SingleElectron_RunH_Prompt_v2,
	SingleElectron_RunH_Prompt_v3,

	SingleMuon_RunB,
	SingleMuon_RunC,
	SingleMuon_RunD,
	SingleMuon_RunE,
	SingleMuon_RunF,
	SingleMuon_RunG,
	SingleMuon_RunH_Prompt_v2,
	SingleMuon_RunH_Prompt_v3,


    TTJets_PowhegPythia8,
    TTJets_PowhegPythia8_fsrdown,
    TTJets_PowhegPythia8_fsrup,
    TTJets_PowhegPythia8_isrdown,
    TTJets_PowhegPythia8_isrup,
    TTJets_PowhegPythia8_down,
    TTJets_PowhegPythia8_up,
    TTJets_PowhegPythia8_mtop1695,
    TTJets_PowhegPythia8_mtop1755,

	TTJets_amcatnloFXFX,
	TTJets_madgraphMLM,
	TTJets_PowhegHerwigpp,

	// TTJets_synch,
	DYJetsToLL_M50,
	DY1JetsToLL_M50,
	DY2JetsToLL_M50,
	DY3JetsToLL_M50,
	DY4JetsToLL_M50,
	WJetsToLNu,
	W1JetsToLNu,
	W2JetsToLNu,
	W3JetsToLNu,
	W4JetsToLNu,
	// TToLeptons_s,
	// TBarToLeptons_s,
	// TBarToLeptons_t,
	ST_t,
	STbar_t,
	ST_s,
	ST_tW,
	STbar_tW,
	// GJets_100to200,
	// GJets_200to400,
	// GJets_400to600,
	// GJets_600toInf,
	QCD_EMEnriched_15to20,
	QCD_EMEnriched_20to30,
	QCD_EMEnriched_30to50,
	QCD_EMEnriched_50to80,
	QCD_EMEnriched_80to120,
	QCD_EMEnriched_120to170,
	QCD_EMEnriched_170to300,
	QCD_EMEnriched_300toInf,
	QCD_bcToE_15to20,
	QCD_bcToE_20to30,
	QCD_bcToE_30to80,
	QCD_bcToE_80to170,
	QCD_bcToE_170to250,
	QCD_bcToE_250toInf,
	QCD_MuEnriched_15to20,
	QCD_MuEnriched_20to30,
	QCD_MuEnriched_30to50,
	QCD_MuEnriched_50to80,
	QCD_MuEnriched_80to120,
	QCD_MuEnriched_120to170,
	QCD_MuEnriched_170to300,
	QCD_MuEnriched_300to470,
	QCD_MuEnriched_470to600,
	QCD_MuEnriched_600to800,
	QCD_MuEnriched_800to1000,
	QCD_MuEnriched_1000toInf,

	QCD_50to80,
	QCD_80to120,
	QCD_120to170,
	QCD_170to300,
	QCD_300to470,
	QCD_470to600,
	QCD_600to800,
	QCD_800to1000,
	// QCD_1000toInf,

	NUMBER_OF_DATA_TYPES
};

const boost::array<std::string, DataType::NUMBER_OF_DATA_TYPES> names = { {
//
		"ElectronHad",//
		// "MuHad", //

		"SingleElectron_RunB",
		"SingleElectron_RunC",
		"SingleElectron_RunD",
		"SingleElectron_RunE",
		"SingleElectron_RunF",
		"SingleElectron_RunG",
		"SingleElectron_RunH_Prompt_v2",
		"SingleElectron_RunH_Prompt_v3",

		"SingleMuon_RunB",
		"SingleMuon_RunC",
		"SingleMuon_RunD",
		"SingleMuon_RunE",
		"SingleMuon_RunF",
		"SingleMuon_RunG",
		"SingleMuon_RunH_Prompt_v2",
		"SingleMuon_RunH_Prompt_v3",

	    "TTJets_PowhegPythia8",
	    "TTJets_PowhegPythia8_fsrdown",
	    "TTJets_PowhegPythia8_fsrup",
	    "TTJets_PowhegPythia8_isrdown",
	    "TTJets_PowhegPythia8_isrup",
	    "TTJets_PowhegPythia8_down",
	    "TTJets_PowhegPythia8_up",
	    "TTJets_PowhegPythia8_mtop1695",
	    "TTJets_PowhegPythia8_mtop1755",

		"TTJets_amcatnloFXFX",
		// "TTJets_amcatnloFXFX_scaledown",
		// "TTJets_amcatnloFXFX_scaleup",
		// "TTJets_amcatnloFXFX_mtop1695",
		// "TTJets_amcatnloFXFX_mtop1755",

		"TTJets_madgraphMLM",
		"TTJets_PowhegHerwigpp",

		// "TTJets_synch",

		"DYJetsToLL_M50",
		"DY1JetsToLL_M50",
		"DY2JetsToLL_M50",
		"DY3JetsToLL_M50",
		"DY4JetsToLL_M50",
		"WJetsToLNu",
		"W1JetsToLNu",
		"W2JetsToLNu",
		"W3JetsToLNu",
		"W4JetsToLNu",
		"ST_t",
		"STbar_t",
		"ST_s",
		"ST_tW",
		"STbar_tW",
		// "GJets_100to200",
		// "GJets_200to400",
		// "GJets_400to600",
		// "GJets_600toInf",
		"QCD_EMEnriched_15to20",
		"QCD_EMEnriched_20to30",
		"QCD_EMEnriched_30to50",
		"QCD_EMEnriched_50to80",
		"QCD_EMEnriched_80to120",
		"QCD_EMEnriched_120to170",
		"QCD_EMEnriched_170to300",
		"QCD_EMEnriched_300toInf",
		"QCD_bcToE_15to20",
		"QCD_bcToE_20to30",
		"QCD_bcToE_30to80",
		"QCD_bcToE_80to170",
		"QCD_bcToE_170to250",
		"QCD_bcToE_250toInf",
		"QCD_MuEnriched_15to20",
		"QCD_MuEnriched_20to30",
		"QCD_MuEnriched_30to50",
		"QCD_MuEnriched_50to80",
		"QCD_MuEnriched_80to120",
		"QCD_MuEnriched_120to170",
		"QCD_MuEnriched_170to300",
		"QCD_MuEnriched_300to470",
		"QCD_MuEnriched_470to600",
		"QCD_MuEnriched_600to800",
		"QCD_MuEnriched_800to1000",
		"QCD_MuEnriched_1000toInf",

		"QCD_50to80",
		"QCD_80to120",
		"QCD_120to170",
		"QCD_170to300",
		"QCD_300to470",
		"QCD_470to600",
		"QCD_600to800",
		"QCD_800to1000",
		// "QCD_1000toInf",
		} };

DataType::value getDataType(const std::string filename);
}

}

#endif /* ENUMERATORS_H_ */
