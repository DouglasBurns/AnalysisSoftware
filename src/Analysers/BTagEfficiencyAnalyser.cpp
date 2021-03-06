/*
 * BTagEff.cpp
 *
 *  Created on: 23 Jul 2015
 *      Author: burns
 */

#include "../../interface/Analysers/BTagEfficiencyAnalyser.h"

using namespace BAT;

void BTagEff::analyse(const EventPtr event) {

	treeMan_->setCurrentFolder(histogramFolder_);
	int NJets = 0;
	int NBJets = 0; // How many medium b jets
	const JetCollection cleanedJets = event->CleanedJets();

	int selectionCriteria = -1;
	if ( event->PassesElectronTriggerAndSelectionNoB() ) selectionCriteria = SelectionCriteria::ElectronPlusJetsReference;
	else if ( event->PassesMuonTriggerAndSelectionNoB() ) selectionCriteria = SelectionCriteria::MuonPlusJetsReference;
	const LeptonPointer signalLepton = event->getSignalLepton( selectionCriteria );

	// unsigned int nParton = 0;
	for (unsigned int jetIndex = 0; jetIndex < cleanedJets.size(); ++jetIndex) {
		const JetPointer jet(cleanedJets.at(jetIndex));

		bool isLoose = false;
		bool isMedium = false;
		bool isTight = false;

		double jetPt = jet->pt();
		double jetEta = jet->eta();

		if (jetPt < 30 || fabs(jetEta) > 2.4) continue;

		double jetCSV = jet->getBTagDiscriminator(BAT::BtagAlgorithm::value::CombinedSecondaryVertexV2);
		if (jet->isLoose()){ 
			isLoose = true;
		}
		if (jet->isMedium()){
			isMedium = true;
			++NBJets;
		} 
		if (jet->isTight()){
			isTight = true;
		}

		unsigned int partonFlavour = abs(jet->partonFlavour());
		unsigned int hadronFlavour = abs(jet->hadronFlavour());

		treeMan_->Fill("pt", jetPt);
		treeMan_->Fill("eta", jetEta);
		treeMan_->Fill("CSV", jetCSV);
		treeMan_->Fill("partonFlavour", partonFlavour);
		treeMan_->Fill("hadronFlavour", hadronFlavour);
		treeMan_->Fill("isLoose", isLoose);
		treeMan_->Fill("isMedium", isMedium);
		treeMan_->Fill("isTight", isTight);
		++NJets;
	}

	treeMan_->Fill("NJets", NJets);
	treeMan_->Fill("NBJets", NBJets);
	treeMan_->Fill("EventWeight", event->weight());
	treeMan_->Fill("PUWeight", event->PileUpWeight());

	if ( selectionCriteria == SelectionCriteria::ElectronPlusJetsReference ) {
		double electronEfficiencyCorrection = 1;
		if ( !event->isRealData() ) {
			const ElectronPointer signalElectron(boost::static_pointer_cast<Electron>(signalLepton));
			electronEfficiencyCorrection = signalElectron->getEfficiencyCorrection( 0 );
		}
		treeMan_->Fill("ElectronEfficiencyCorrection",electronEfficiencyCorrection);
	}
	else if ( selectionCriteria == SelectionCriteria::MuonPlusJetsReference ) {
		double muonEfficiencyCorrection = 1;
		if ( !event->isRealData() ) {
			const MuonPointer signalMuon(boost::static_pointer_cast<Muon>(signalLepton));
			muonEfficiencyCorrection = signalMuon->getEfficiencyCorrection( 0 );
		}
		treeMan_->Fill("MuonEfficiencyCorrection",muonEfficiencyCorrection);
	}
}


BTagEff::BTagEff(TreeManagerPtr treeMan, std::string histogramFolder) :
	BasicAnalyser(treeMan, histogramFolder){
}

BTagEff::~BTagEff() {
}

void BTagEff::createTrees() {
	treeMan_->setCurrentFolder(histogramFolder_);

	treeMan_->addBranch("pt", "F", "Jets" + Globals::treePrefix_, false);
	treeMan_->addBranch("eta", "F", "Jets" + Globals::treePrefix_, false);
	treeMan_->addBranch("CSV", "F", "Jets" + Globals::treePrefix_, false);
	treeMan_->addBranch("partonFlavour", "F", "Jets" + Globals::treePrefix_, false);
	treeMan_->addBranch("hadronFlavour", "F", "Jets" + Globals::treePrefix_, false);
	treeMan_->addBranch("isLoose", "F", "Jets" + Globals::treePrefix_, false);
	treeMan_->addBranch("isMedium", "F", "Jets" + Globals::treePrefix_, false);
	treeMan_->addBranch("isTight", "F", "Jets" + Globals::treePrefix_, false);
	treeMan_->addBranch("NJets", "F", "Jets" + Globals::treePrefix_);
	treeMan_->addBranch("NBJets", "F", "Jets" + Globals::treePrefix_);
	treeMan_->addBranch("ElectronEfficiencyCorrection", "F", "Jets" + Globals::treePrefix_);
	treeMan_->addBranch("MuonEfficiencyCorrection", "F", "Jets" + Globals::treePrefix_);
}
/* namespace BAT */
