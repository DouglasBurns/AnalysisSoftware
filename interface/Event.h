/*
 * Event.h
 *
 *  Created on: Jun 25, 2010
 *      Author: lkreczko
 */

#ifndef EVENT_H_
#define EVENT_H_
#include "RecoObjects/Electron.h"
#include "RecoObjects/Jet.h"
#include "RecoObjects/Muon.h"
#include "RecoObjects/Vertex.h"
#include "RecoObjects/MET.h"
#include "RecoObjects/MCParticle.h"
#include "RecoObjects/PseudoTopParticles.h"
#include "RecoObjects/SelectionOutputInfo.h"
#include "RecoObjects/TTGenInfo.h"
#include "DataTypes.h"
#include "Printers/EventContentPrinter.h"
#include "HighLevelTriggers.h"
#include "TtbarHypothesis.h"

#include <boost/shared_ptr.hpp>
#include <string>
#include <bitset>
 
//forward definition of event needed before BasicSelection inclusion.
namespace BAT {
class Event;
typedef boost::shared_ptr<Event> EventPtr;
}
#include "GlobalVariables.h"
namespace BAT {

namespace SelectionCriteria {
enum selection {
	MuonPlusJetsReference = 1,
	ElectronPlusJetsReference = 2,
	MuonPlusJetsQCDNonIsolated1p5to3 = 3,
	MuonPlusJetsQCDNonIsolated3toInf = 4,
	ElectronPlusJetsQCDNonIsolated = 5,
	ElectronPlusJetsQCDConversion = 6,
	NUMBER_OF_SELECTION_STEPS
};
}

class Event {
protected:
	boost::shared_ptr<std::vector<int> > HLTs;
	boost::shared_ptr<std::vector<int> > HLTPrescales;
	VertexCollection vertices;
	VertexCollection goodVertices;
	TrackCollection tracks;
	ElectronCollection allElectrons;

	JetCollection allJets;
	JetCollection cleanedJets_;
	JetCollection cleanedBJets_;
	JetCollection genJets;

	MuonCollection allMuons;

	bool passesElectronChannelTrigger_;
	bool passesElectronChannelQCDTrigger_;
	bool passesMuonChannelTrigger_;
	bool passesTkMuonChannelTrigger_;
	bool passesMuonChannelQCDTrigger_;

	bool passesElectronSelection_;
	bool passesElectronQCDSelection_;
	bool passesElectronConversionSelection_;
	bool passesMuonSelection_;
	bool passesMuonQCDSelection1p5to3_;
	bool passesMuonQCDSelection3toInf_;

	bool passesElectronSelectionNoB_;
	bool passesElectronQCDSelectionNoB_;
	bool passesElectronConversionSelectionNoB_;
	bool passesMuonSelectionNoB_;
	bool passesMuonQCDSelection1p5to3NoB_;
	bool passesMuonQCDSelection3toInfNoB_;

	bool passesMETFilters_;

	bool isSemiLeptonicElectron_;
	bool isSemiLeptonicMuon_;
	SelectionOutputInfo selectionOutputInfo_electron;
	SelectionOutputInfo selectionOutputInfo_muon;
	SelectionOutputInfo selectionOutputInfo_electronQCDNonisolated;
	SelectionOutputInfo selectionOutputInfo_electronQCDConversion;
	SelectionOutputInfo selectionOutputInfo_muonQCDNonisolated1p5to3;
	SelectionOutputInfo selectionOutputInfo_muonQCDNonisolated3toInf;

	TTGenInfoPointer ttGenInfo_;

	PseudoTopParticlesPointer pseudoTopParticles_;
	MCParticleCollection genParticles;

	METCollection mets_;
	METPointer genMet_;

	TtbarHypothesis ttbarHypothesis_;

	DataType::value dataType;
	unsigned long runNumber;
	unsigned long eventNumber;
	unsigned long localEventNumber;
	unsigned long lumiBlock;
	double eventWeight;
	double pileUpWeight_, pileUpWeight_up_, pileUpWeight_down_;
	double bJetWeight_, lightJetUpWeight_, lightJetDownWeight_, bJetUpWeight_, bJetDownWeight_;
	double bJetEfficiencyCorrectionWeight_, bJetAlternativeWeight_;
	double tightbJetWeight_, tightlightJetUpWeight_, tightlightJetDownWeight_, tightbJetUpWeight_, tightbJetDownWeight_;
	//These values are calculated with PileUp_2011_truth_finebin.root
	double PUWeightInTimeOnly_, PUWeight3BX_, // PUWeight3D_,
			PUWeightShiftUp_, PUWeightShiftDown_;

	unsigned int numberOfHighPurityTracks_;
	bool isBeamScraping_;

	std::vector<int> genNumberOfPileUpVertices;
	std::vector<int> trueNumberOfPileUpVertices_;
	double generatorWeight_;
	double centralLHEWeight_;
	std::vector<double> generatorSystematicWeights_;
	std::vector<int> generatorSystematicWeightsID_;
	double semiLepBrUpWeight_;
	double semiLepBrDownWeight_;
	double petersonFragWeight_;
	double centralFragWeight_;
	double upFragWeight_;
	double downFragWeight_;
	float tau1_, tau2_, tau3_, tau4_;

	double ptdensityRho;
	std::string file_;

	bool passesCSCTightBeamHaloFilter_, passesHBHENoiseFilter_, passesHCALLaserFilter_, passesECALDeadCellFilter_, passesECALDeadCellTPFilter_,
			passesTrackingFailureFilter_, passesEEBadSCFilter_, passesECALLaserCorrFilter_, passesTrackingPOGFilters_;

public:
	Event();
	virtual ~Event();
	bool isRealData() const;
	const DataType::value getDataType() const;
	bool isTTJet(DataType::value type) const;
	void setDataType(DataType::value type);
	void setVertices(VertexCollection vertices);
	void setTracks(TrackCollection tracks);
	void setGenParticles(MCParticleCollection genParticles);
	void setElectrons(ElectronCollection electrons);
	void setJets(JetCollection jets);
	void setJetTTBarPartons();
	void setCleanedJets(JetCollection jets);
	void setCleanedBJets(JetCollection bjets);
	void setGenJets(JetCollection genJets);
	void setMuons(MuonCollection muons);
	void setPassesElectronChannelTrigger( bool passesTrigger );
	void setPassesElectronChannelQCDTrigger( bool passesTrigger );
	void setPassesMuonChannelTrigger( bool passesTrigger );
	void setPassesTkMuonChannelTrigger( bool passesTrigger );
	void setPassesMuonChannelQCDTrigger( bool passesTrigger );
	void setPassesElectronSelection(bool passesElectronSelection);
	void setPassesElectronQCDSelection(bool passesElectronQCDSelection);
	void setPassesElectronConversionSelection(bool passesElectronConversionSelection);
	void setPassesMuonSelection(bool passesMuonSelection);
	void setPassesMuonQCDSelection1p5to3(bool passesMuonQCDSelection1p5to3);
	void setPassesMuonQCDSelection3toInf(bool passesMuonQCDSelection3toInf);
	void setPassesElectronSelectionNoB(bool passesElectronSelectionNoB);
	void setPassesElectronQCDSelectionNoB(bool passesElectronQCDSelectionNoB);
	void setPassesElectronConversionSelectionNoB(bool passesElectronConversionSelectionNoB);
	void setPassesMuonSelectionNoB(bool passesMuonSelectionNoB);
	void setPassesMuonQCDSelection1p5to3NoB(bool passesMuonQCDSelection1p5to3NoB);
	void setPassesMuonQCDSelection3toInfNoB(bool passesMuonQCDSelection3toInfNoB);
	void setPassOfflineSelectionInfo( std::vector<unsigned int> );
	const bool passesJetSelection( const unsigned int selectionCriteria );
	const bool passesBJetSelection( const unsigned int selectionCriteria );
	const bool passesSignalLeptonSelection( const unsigned int selectionCriteria );
	void setPassGenSelectionInfo( std::vector<unsigned int> );
	void setIsSemiLeptonicElectron( bool isSemiLeptonicElectron );
	void setIsSemiLeptonicMuon( bool isSemiLeptonicMuon );
	void setElectronSelectionOutputInfo(SelectionOutputInfo newSelectionOutputInfo);
	void setMuonSelectionOutputInfo(SelectionOutputInfo newSelectionOutputInfo);
	void setElectronQCDNonisolatedSelectionOutputInfo(SelectionOutputInfo newSelectionOutputInfo);
	void setElectronConversionSelectionOutputInfo(SelectionOutputInfo newSelectionOutputInfo);
	void setMuonQCDNonisolatedSelection1p5to3OutputInfo(SelectionOutputInfo newSelectionOutputInfo);
	void setMuonQCDNonisolatedSelection3toInfOutputInfo(SelectionOutputInfo newSelectionOutputInfo);
	void setTTGenInfo(TTGenInfoPointer ttGenInfo );
	void setPseudoTopParticles(PseudoTopParticlesPointer newPseudoParticles );
	void setMETs(const std::vector<METPointer> mets);
	void setGenMET(const METPointer met);
	void setTTbarHypothesis(const TtbarHypothesis newHypo);
	void setHLTs(const boost::shared_ptr<std::vector<int> >);
	void setHLTPrescales(const boost::shared_ptr<std::vector<int> >);
	void setFile(std::string file);
	void setRunNumber(unsigned long number);
	void setEventNumber(unsigned long number);
	void setLocalEventNumber(unsigned long number);
	void setPassesMETFilters( bool pass );
	void setLumiBlock(unsigned long block);
	void setEventWeight(double weight);
	void setPileUpWeight(double weight, int systematic=0);
	void setBJetWeight(double weight);
	void setBJetEfficiencyCorrectionWeight(double weight);
	void setBJetAlternativeWeight(double weight);
	void setLightJetUpWeight(double weight);
	void setLightJetDownWeight(double weight);
	void setBJetUpWeight(double weight);
	void setBJetDownWeight(double weight);
	void setTightBJetWeight(double weight);
	void setTightLightJetUpWeight(double weight);
	void setTightLightJetDownWeight(double weight);
	void setTightBJetUpWeight(double weight);
	void setTightBJetDownWeight(double weight);
	void setBeamScrapingVeto(bool isScraping);
	void setGenNumberOfPileUpVertices(std::vector<int> pileup);
	void setTrueNumberOfPileUpVertices(std::vector<int> pileup);
	void setGeneratorWeight( double generatorWeight );
	void setCentralLHEWeight( double centralLHEWeight );
	void setGeneratorSystematicWeights( std::vector<double> generatorSystematicWeights );
	void setGeneratorSystematicWeightsID( std::vector<int> generatorSystematicWeightsID );
	void setSemiLepBrUpWeight( double weight );
	void setSemiLepBrDownWeight( double weight );
	void setPetersonFragWeight( double weight );
	void setCentralFragWeight( double weight );
	void setUpFragWeight( double weight );
	void setDownFragWeight( double weight );
	void setTau1( float tau );
	void setTau2( float tau );
	void setTau3( float tau );
	void setTau4( float tau );
	void setPtDensityRho(double rho);
	void setPUWeightInTimeOnly(double puweight);
	void setPUWeight3BX(double puweight);
	void setPUWeightShiftUp(double puweight);
	void setPUWeightShiftDown(double puweight);
	void setCSCTightBeamHaloFilter(bool result);
	void setHBHENoiseFilter(bool result);
	void setHCALLaserFilter(bool result);
	void setECALDeadCellFilter(bool result);
	void setECALDeadCellTPFilter(bool result);
	void setTrackingFailureFilter(bool result);
	void setEEBadSCFilter(bool result);
	void setECALLaserCorrFilter(bool result);
	void setTrackingPOGFilters(bool result);

	const VertexPointer PrimaryVertex() const;
	const VertexCollection& Vertices() const;
	const VertexCollection& GoodVertices() const;
	const TrackCollection& Tracks() const;
	const MCParticleCollection& GenParticles() const;
	const ElectronCollection& Electrons() const;
	const ElectronCollection& QCDElectrons() const;
	const JetCollection& Jets() const;
	const JetCollection& CleanedJets() const;
	const JetCollection& CleanedBJets() const;
	const JetCollection& GenJets() const;
	const MuonCollection& Muons() const;
	const METPointer MET() const;
	const METPointer MET(METAlgorithm::value type) const;
	const METPointer GenMET() const;

	const TtbarHypothesis ttbarHypothesis() const;

	const bool PassesElectronChannelTrigger() const;
	const bool PassesMuonChannelTrigger() const;
	const bool PassesTkMuonChannelTrigger() const;

	const bool PassesElectronSelection() const;
	const bool PassesElectronQCDSelection() const;
	const bool PassesElectronConversionSelection() const;
	const bool PassesMuonSelection() const;
	const bool PassesMuonQCDSelection1p5to3() const;
	const bool PassesMuonQCDSelection3toInf() const;

	const bool PassesElectronTriggerAndSelection() const;
	const bool PassesMuonTriggerAndSelection() const;
	const bool PassesElectronTriggerAndQCDSelection() const;
	const bool PassesElectronTriggerAndConversionSelection() const;
	const bool PassesMuonTriggerAndQCDSelection1p5to3() const;
	const bool PassesMuonTriggerAndQCDSelection3toInf() const;
		
	const bool PassesElectronSelectionNoB() const;
	const bool PassesElectronQCDSelectionNoB() const;
	const bool PassesElectronConversionSelectionNoB() const;
	const bool PassesMuonSelectionNoB() const;
	const bool PassesMuonQCDSelection1p5to3NoB() const;
	const bool PassesMuonQCDSelection3toInfNoB() const;

	const bool PassesElectronTriggerAndSelectionNoB() const;
	const bool PassesMuonTriggerAndSelectionNoB() const;
	const bool PassesElectronTriggerAndQCDSelectionNoB() const;
	const bool PassesElectronTriggerAndConversionSelectionNoB() const;
	const bool PassesMuonTriggerAndQCDSelection1p5to3NoB() const;
	const bool PassesMuonTriggerAndQCDSelection3toInfNoB() const;

	const bool isSemiLeptonicElectron() const;
	const bool isSemiLeptonicMuon() const;
	
	const TTGenInfoPointer TTGenInfo() const;

	const PseudoTopParticlesPointer PseudoTopParticles() const;

	const LeptonPointer getSignalLepton( unsigned int selectionCriteria ) const;
	const JetCollection getCleanedJets( unsigned int selectionCriteria ) const;
	const JetCollection getCleanedBJets( unsigned int selectionCriteria ) const;
	const unsigned int getNBJets( unsigned int selectionCriteria ) const;
	const std::vector<bool> getCuts( unsigned int selectionCriteria ) const;

	JetCollection GetBJetCollection(const JetCollection& jets, BtagAlgorithm::value btagAlgorithm,
			BtagAlgorithm::workingPoint WP = BtagAlgorithm::MEDIUM) const;
	std::string file() const;
	unsigned long runnumber() const;
	unsigned long eventnumber() const;
	unsigned long localnumber() const;
	unsigned long lumiblock() const;
	double weight() const;
	double PileUpWeight( int systematic = 0 ) const;

	double BJetWeight() const;
	double BJetEfficiencyCorrectionWeight() const;
	double BJetAlternativeWeight() const;

	double LightJetUpWeight() const;
	double LightJetDownWeight() const;
	double BJetUpWeight() const;
	double BJetDownWeight() const;

	double TightBJetWeight() const;
	double TightLightJetUpWeight() const;
	double TightLightJetDownWeight() const;
	double TightBJetUpWeight() const;
	double TightBJetDownWeight() const;

	void inspect() const;
	bool HLT(HLTriggers::value trigger) const;
	int HLTPrescale(HLTriggers::value trigger) const;
	static bool useCustomConversionTagger;
	static bool usePFIsolation;

	const std::vector<int> GeneratedPileUpVertices() const;
	const double generatorWeight() const;
	const double centralLHEWeight() const;
	const std::vector<double> generatorSystematicWeights() const;
	const std::vector<int> generatorSystematicWeightsID() const;
	const double semiLepBrUpWeight() const;
	const double semiLepBrDownWeight() const;
	const double petersonFragWeight() const;
	const double centralFragWeight() const;
	const double upFragWeight() const;
	const double downFragWeight() const;

	const float tau1() const;
	const float tau2() const;
	const float tau3() const;
	const float tau4() const;

	double averageNumberOfVertices() const;
	double inTimeOnlyNumberOfVertices() const;
	double rho() const;
	unsigned int numberOfHighPurityTracks() const;
	bool isBeamScraping() const;
	const std::vector<int>& getTrueNumberOfVertices() const;

	double PUWeightInTimeOnly() const;
	double PUWeight3BX() const;
	double PUWeightShiftUp() const;
	double PUWeightShiftDown() const;

	bool passesCSCTightBeamHaloFilter() const;
	bool passesHBHENoiseFilter() const;
	bool passesHCALLaserFilter() const;
	bool passesECALDeadCellFilter() const;
	bool passesECALDeadCellTPFilter() const;
	bool passesTrackingFailureFilter() const;
	bool passesEEBadSCFilter() const;
	bool passesECALLaserCorrFilter() const;
	bool passesTrackingPOGFilters() const;

	bool passesMETFilters() const;

	static double HT(const JetCollection jets, const double jetPtThreshold = minJetPt_ );
	static double M3(const JetCollection jets, const double jetPtThreshold = minJetPt_ );
	static double M_bl(const JetCollection, const ParticlePointer, const double jetPtThreshold = minJetPt_ );
	static double angle_bl(const JetCollection, const ParticlePointer, const double jetPtThreshold = minJetPt_ );
	static double pseudo_angle_bl( const MCParticleCollection pseudoBs, const ParticlePointer pseudoLepton, const double jetPtThreshold = minJetPt_ );
	static void getTopTwoCSVJets(const JetCollection, unsigned int& highestCSVJetIndex, unsigned int& secondHighestCSVJetIndex);
	static double deltaPhi_bb(const ParticlePointer, const ParticlePointer);
	static double deltaEta_bb(const ParticlePointer, const ParticlePointer);
	static double angle_bb(const ParticlePointer, const ParticlePointer);
	static unsigned int NJets(const JetCollection, const double jetPtThreshold = minJetPt_ );
	static double ST(const JetCollection, const ParticlePointer, const METPointer, const double jetPtThreshold = minJetPt_ );
	static double MT(const ParticlePointer, const METPointer);
	static double WPT(const ParticlePointer, const METPointer);

	static double const minJetPt_;
	static double const maxJetAbsEta_;
	static unsigned int const minNJets_;
	static unsigned int const minNBJets_;

	static double const minSignalMuonPt_;
	static double const minSignalMuonEta_;
	static double const minSignalElectronPt_;
	static double const minSignalElectronEta_;

};

}

#endif /* EVENT_H_ */
