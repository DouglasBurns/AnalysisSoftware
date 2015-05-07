/*
 * LikelihoodInputAnalyser.cpp
 *
 *  Created on: May 5, 2015
 *      Author: ejclemen
 */
#include "../../interface/Analysers/LikelihoodInputAnalyser.h"
#include "../../interface/GlobalVariables.h"
#include <boost/scoped_ptr.hpp>

using namespace BAT;


void LikelihoodInputAnalyser::analyse(const EventPtr event) {

	weight_ = event->weight();
	treeMan_->setCurrentFolder(histogramFolder_);

	// Require one of the event selections to have been satisfied, and we have a genuine ttbar event
	int selectionCriteria = -1;
	if ( event->PassesElectronTriggerAndSelection() && event->isSemiLeptonicElectron() ) selectionCriteria = SelectionCriteria::ElectronPlusJetsReference;
	else if ( event->PassesMuonTriggerAndSelection() && event->isSemiLeptonicMuon() ) selectionCriteria = SelectionCriteria::MuonPlusJetsReference;

	if ( selectionCriteria < 0 ) return;

	// Get objects
	const JetCollection jets(event->CleanedJets());
	const JetCollection bjets(event->CleanedBJets());
	const LeptonPointer signalLepton = event->getSignalLepton( selectionCriteria );
	METAlgorithm::value metType = (METAlgorithm::value) 0;
	const METPointer met(event->MET(metType));
	const JetCollection genJets( event->GenJets() );

	// Loop over jets and see if the parton matched to it is from ttbar decay
	for ( unsigned int jet1Index=0; jet1Index < jets.size(); ++jet1Index ) {
		JetPointer jet1 = jets[jet1Index];
		cout << "Jet pt : " << jet1->pt() << endl;
		cout << "TTbar parton : " << jet1->ttbar_decay_parton() << endl;
	}

}

LikelihoodInputAnalyser::LikelihoodInputAnalyser(HistogramManagerPtr histMan, boost::shared_ptr<TreeManager> treeMan, std::string histogramFolder ) :
		BasicAnalyser(histMan, treeMan, histogramFolder)
{

}

LikelihoodInputAnalyser::~LikelihoodInputAnalyser() {
}

void LikelihoodInputAnalyser::createHistograms() {
	histMan_->setCurrentHistogramFolder(histogramFolder_);
}

void LikelihoodInputAnalyser::createTrees() {
	treeMan_->setCurrentFolder(histogramFolder_);

	treeMan_->addBranch("BRANCH", "F", "Likelihood" + Globals::treePrefix_);
}
