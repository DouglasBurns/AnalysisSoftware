/*
 * LikelihoodInputAnalyser.cpp
 *
 *  Created on: May 5, 2015
 *      Author: ejclemen
 */
#include "../../interface/Analysers/LikelihoodInputAnalyser.h"
#include "../../interface/GlobalVariables.h"
#include <boost/scoped_ptr.hpp>
#include "../Rochester/NeutrinoSolver.cc" // WHY???

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

	cout << endl << "*****************************************" << endl;
	cout << "Number of Jets : " << jets.size() << endl;
	cout << "Number of bJets : " << bjets.size() << endl;
	cout << "Number of Light Jets : " << jets.size() - bjets.size() << endl;

	// HADRONIC W

	for ( unsigned int jet1Index=0; jet1Index < (jets.size()-1); ++jet1Index ) {
		cout << "------------------" << endl;
		FourVector W;
		// double MW=0;
		JetPointer jet1 = jets[jet1Index];
		if (jet1->ttbar_decay_parton() >= 5) continue;

		for ( unsigned int jet2Index=1; jet2Index < jets.size(); ++jet2Index ) {
			if (jet1Index >= jet2Index) continue;

			JetPointer jet2 = jets[jet2Index];
			if (jet2->ttbar_decay_parton() >= 5) continue;

			W = jet1->getFourVector() + jet2->getFourVector();//fourvector addition = A + B
			// MW = W.M();
			cout << "Mass of W : " << W.M() << "GeV" << endl;
		}
	}

	// // Loop over jets and see if the parton matched to it is from ttbar decay
	// for ( unsigned int jet1Index=0; jet1Index < jets.size(); ++jet1Index ) {
	// 	JetPointer jet1 = jets[jet1Index];
	// 	cout << "Jet pt : " << jet1->pt() << endl;
	// 	cout << "TTbar parton : " << jet1->ttbar_decay_parton() << endl; 
	// }

	// Loop over b jets calculate is Good/Bad and store chi sq neutrino
	// cout << "-----------------------------------------" << endl;

	for ( unsigned int bjetIndex=0; bjetIndex < bjets.size(); ++bjetIndex ) {

		JetPointer bjet = bjets[bjetIndex];
		bool leptonicB = false;
		bool hadronicB = false;

		// cout << "-----------------------------------------" << endl;
		// cout << "Bjet parton : " << bjet->ttbar_decay_parton() << endl;//to see whether b jet

		if (bjet->ttbar_decay_parton() == 5)
			leptonicB = true;
		if (bjet->ttbar_decay_parton() == 6)
			hadronicB = true;

		NeutrinoSolver neutrinoSolver( &(signalLepton->getFourVector()), &(bjet->getFourVector()), 80, 173 );
		double test = -1;
		FourVector neutrino = neutrinoSolver.GetBest(met->getFourVector().Px(), met->getFourVector().Py(), 25., 25., 0., test );
		// cout << "Neutrino ChiSq : " << test << endl;

		if (leptonicB == true)
			treeMan_->Fill("Leptonic_B", test);
		else
			treeMan_->Fill("False_B", test);
			if  (hadronicB == true)
				treeMan_->Fill("Hadronic_B", test);

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
	treeMan_->addBranch("Leptonic_B", "F", "Neutrino" + Globals::treePrefix_);
	treeMan_->addBranch("Hadronic_B", "F", "Neutrino" + Globals::treePrefix_);
	treeMan_->addBranch("False_B", "F", "Neutrino" + Globals::treePrefix_);

}
