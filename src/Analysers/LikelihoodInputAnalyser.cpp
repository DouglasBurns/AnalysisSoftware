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
	unsigned int parton = 0;
	double CSV = -1;
	int TypeOfJet = 0;

	const JetCollection AllJets(event->Jets());
	const JetCollection jets(event->CleanedJets());
	const JetCollection bjets(event->CleanedBJets());

	METAlgorithm::value metType = (METAlgorithm::value) 0;
	const METPointer met(event->MET(metType));

	if (event->isSemiLeptonicElectron() || event->isSemiLeptonicMuon()){

		for ( unsigned int AllJetIndex=0; AllJetIndex < AllJets.size(); ++AllJetIndex ) {

			JetPointer AllJet = AllJets[AllJetIndex];

			if (AllJet->matched_parton() == 0) continue;//if null pointer

			parton = abs(AllJet->matched_parton()->pdgId());
			CSV = AllJet->getBTagDiscriminator(BAT::BtagAlgorithm::value::CombinedSecondaryVertexV2);

			if ( (parton >= 1 && parton <= 4) || (parton == 21) ){
				TypeOfJet = 1;
				treeMan_->Fill("TypeofJet", TypeOfJet);
			}

			else if ( parton == 5 ){
				TypeOfJet = 2;
				treeMan_->Fill("TypeofJet", TypeOfJet);
			}
			treeMan_->Fill("Jets", CSV);
		}
	}

	// Require one of the event selections to have been satisfied, and we have a genuine ttbar event
	int selectionCriteria = -1;
	if ( event->PassesElectronTriggerAndSelection() && event->isSemiLeptonicElectron() ) selectionCriteria = SelectionCriteria::ElectronPlusJetsReference;
	else if ( event->PassesMuonTriggerAndSelection() && event->isSemiLeptonicMuon() ) selectionCriteria = SelectionCriteria::MuonPlusJetsReference;
	// cout << "Pass Electron Trig & Sel : " << event->PassesElectronTriggerAndSelection() << " Pass Muon Trig & Sel : " << event->PassesMuonTriggerAndSelection() << endl;
	// cout << "Pass is emileptonic electron : " << event->isSemiLeptonicElectron() << " Pass semileptonic muon : " << event->isSemiLeptonicMuon() << endl;
	// cout << "selection criteria : " << selectionCriteria << endl;
	if ( selectionCriteria < 0 ) return;
	// cout << "***********************" << endl;
	// cout << "All Jets Size : " << AllJets.size() << endl;
	// cout << "B Jets Size : " << bjets.size() << endl;
	// cout << "Jets Size : " << jets.size() << endl;

	// Get objects
	const LeptonPointer signalLepton = event->getSignalLepton( selectionCriteria );

	JetCollection jetsWithoutBs;
	FourVector WLep, TopLep, WHad, TopHad;

	// Get cleaned jets that aren't b tagged
	for ( unsigned int jetIndex=0; jetIndex < jets.size(); ++jetIndex ) {
		bool isBJet = false;
		JetPointer thisJet = jets[jetIndex];
		for ( unsigned int bJetIndex=0; bJetIndex < bjets.size(); ++bJetIndex ) {
			JetPointer thisBJet = bjets[bJetIndex];
			if ( thisJet == thisBJet ) {
				isBJet = true;
				break;
			}
		}
		if ( !isBJet ) jetsWithoutBs.push_back( thisJet );
	}

	// cout << "Jets Without B : " << jetsWithoutBs.size() << endl;

	if ( jetsWithoutBs.size() < 2 ) return;

	// cout << endl << "*****************************************" << endl;
	// cout << "Number of Jets : " << jets.size() << endl;
	// cout << "Number of bJets : " << bjets.size() << endl;
	// cout << "Number of Light Jets : " << jets.size() - bjets.size() << endl;
	// cout << "-----------------------------------------" << endl;
	bool Firsttime = true;

	for ( unsigned int bjetIndex=0; bjetIndex < bjets.size(); ++bjetIndex ) {

		JetPointer bjet = bjets[bjetIndex];
		bool leptonicB = false;
		bool mistaggedB = false;
		double CorrectEvent = 0;
		double RecoNu = 0;

		if ( bjet->matched_parton() == 0 ) continue;//if null pointer

		if (bjet->ttbar_decay_parton() == 5){
			leptonicB = true;
			// cout << "B is Actually Leptonic: " << endl;
		}

		if (bjet->ttbar_decay_parton() != 5){
			mistaggedB = true;
			// cout << "B is Mis-Tagged: " << endl;
		}

		if (Firsttime == true){
			for ( unsigned int jetIndex=0; jetIndex < jetsWithoutBs.size(); ++jetIndex ) {
				JetPointer jet = jetsWithoutBs[jetIndex];
				RecoNu = 3;//neutrino comes from normal jet.

				NeutrinoSolver neutrinoSolver( &(signalLepton->getFourVector()), &(jet->getFourVector()), 80, 173 );
				double test = -1;
				FourVector neutrino = neutrinoSolver.GetBest(met->getFourVector().Px(), met->getFourVector().Py(), 25., 25., 0., test );

				treeMan_->Fill("NuChiSq", test);
				treeMan_->Fill("ReconstructedNeutrino", RecoNu);
				// cout << "W " << "Neutrino ChiSq : " << test << endl;
			}
		}
		Firsttime = false;


		//For particular b jet and signal lepton - return the best solution for the neutrino given the t/W mass constraints
		NeutrinoSolver neutrinoSolver( &(signalLepton->getFourVector()), &(bjet->getFourVector()), 80, 173 );
		double test = -1;
		FourVector neutrino = neutrinoSolver.GetBest(met->getFourVector().Px(), met->getFourVector().Py(), 25., 25., 0., test );
		// cout << "Neutrino ChiSq : " << test << endl;
		treeMan_->Fill("NuChiSq", test);

		if (leptonicB == true){
			RecoNu = 1;
			treeMan_->Fill("ReconstructedNeutrino", RecoNu);
		}

		if (mistaggedB == true){
			RecoNu = 2;
			treeMan_->Fill("ReconstructedNeutrino", RecoNu);
		}
		// Calculate the leptonic W and hence leptonic top four vectors.
		WLep = neutrino + signalLepton->getFourVector();
		TopLep = WLep + bjet->getFourVector();

		// cout << "-----------------------------------------" << endl;
		// cout << "Mass of Leptonic W : " << WLep.M() << "GeV" << endl;
		// cout << "Mass of Leptonic Top : " << TopLep.M() << "GeV" << endl;

		//Store the outcome of the Chisq test for leptonic/hadronic and incorrect b jets (inc hadronic)



		// Reconstruct the hadronic W and hence top.
		for ( unsigned int jet1Index=0; jet1Index < (jetsWithoutBs.size()-1); ++jet1Index ) {
			for ( unsigned int jet2Index=jet1Index+1; jet2Index < jetsWithoutBs.size(); ++jet2Index ) {
				JetPointer jet1 = jetsWithoutBs[jet1Index];
				JetPointer jet2 = jetsWithoutBs[jet2Index];

				WHad = jet1->getFourVector() + jet2->getFourVector();
				// cout << "-----------------------------------------" << endl;
				// cout << "Mass of Hadronic W : " << WHad.M() << "GeV" << endl;

				for ( unsigned int AltbjetIndex=0; AltbjetIndex < bjets.size(); ++AltbjetIndex ) {
					if (AltbjetIndex == bjetIndex) continue;
					JetPointer AltBjet = bjets[AltbjetIndex];

					TopHad = WHad + AltBjet->getFourVector();
					// cout << "Mass of Hadronic Top : " << TopHad.M() << "GeV" << endl;

					if ((jet1->ttbar_decay_parton() == 3 && jet2->ttbar_decay_parton() == 4) || (jet2->ttbar_decay_parton() == 3 && jet1->ttbar_decay_parton() == 4)){
						if (AltBjet->ttbar_decay_parton() == 6){
							CorrectEvent = 1;
							if (leptonicB == true){
								CorrectEvent = 2;
							}
						}
					}

					treeMan_->Fill("LeptonicTop_Mass",TopLep.M() );
					treeMan_->Fill("LeptonicW_Mass",WLep.M() );
					treeMan_->Fill("HadronicTop_Mass",TopHad.M() );
					treeMan_->Fill("HadronicW_Mass",WHad.M() );
					treeMan_->Fill("CorrectReconstruction", CorrectEvent);
				}
			}
		}
	}
}

LikelihoodInputAnalyser::LikelihoodInputAnalyser(HistogramManagerPtr histMan, boost::shared_ptr<TreeManager> treeMan, std::string histogramFolder ) :
		BasicAnalyser(histMan, treeMan, histogramFolder){
}

LikelihoodInputAnalyser::~LikelihoodInputAnalyser() {
}

void LikelihoodInputAnalyser::createHistograms() {
	histMan_->setCurrentHistogramFolder(histogramFolder_);
}

void LikelihoodInputAnalyser::createTrees() {
	treeMan_->setCurrentFolder(histogramFolder_);

	treeMan_->addBranch("BRANCH", "F", "Likelihood" + Globals::treePrefix_);

	treeMan_->addBranch("NuChiSq", "F", "Neutrino" + Globals::treePrefix_, false);//false stores vector of event i.e. nuchi2 of all jets - not just last one.
	treeMan_->addBranch("ReconstructedNeutrino", "F", "Neutrino" + Globals::treePrefix_, false);

	treeMan_->addBranch("Jets", "F", "CSV" + Globals::treePrefix_, false);
	treeMan_->addBranch("TypeofJet", "F", "CSV" + Globals::treePrefix_,false);

	treeMan_->addBranch("LeptonicTop_Mass", "F", "TopReco" + Globals::treePrefix_, false);
	treeMan_->addBranch("HadronicTop_Mass", "F", "TopReco" + Globals::treePrefix_, false);
	treeMan_->addBranch("LeptonicW_Mass", "F", "TopReco" + Globals::treePrefix_, false);
	treeMan_->addBranch("HadronicW_Mass", "F", "TopReco" + Globals::treePrefix_, false);
	treeMan_->addBranch("CorrectReconstruction", "F", "TopReco"+ Globals::treePrefix_, false);
}
