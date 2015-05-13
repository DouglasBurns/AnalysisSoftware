/*
 * LikelihoodRecoAnalyser.cpp
 *
 *  Created on: Mey 13, 2015
 *      Author: ejclemen
 */
#include "../../interface/Analysers/LikelihoodRecoAnalyser.h"
#include "../../interface/LumiReWeighting.h"
#include "../../interface/EventWeightProvider.h"
#include "../../interface/PoissonMeanShifter.h"
#include "../../interface/GlobalVariables.h"
#include <boost/scoped_ptr.hpp>

//using namespace reweight;
using namespace BAT;

void LikelihoodRecoAnalyser::analyse(const EventPtr event) {
}

void LikelihoodRecoAnalyser::analyse(const EventPtr event, const JetCollection jets, const JetCollection bjets, const LeptonPointer selectedLepton, const METPointer met ) {

	weight_ = event->weight();
	treeMan_->setCurrentFolder(histogramFolder_);
	histMan_->setCurrentHistogramFolder(histogramFolder_);

	LikelihoodBasedTopPairReconstruction likelihoodReco( selectedLepton, met, jets, bjets );

	likelihoodReco.reconstruct();

	TtbarHypothesisCollection allSolutions = likelihoodReco.getAllSolutions();
	cout << "Got " << allSolutions.size() << " solutions" << endl;

}

LikelihoodRecoAnalyser::LikelihoodRecoAnalyser(HistogramManagerPtr histMan, TreeManagerPtr treeMan, const bool isElectronChannel, std::string histogramFolder ) :
		BasicAnalyser(histMan, treeMan, histogramFolder) //
{

}

LikelihoodRecoAnalyser::~LikelihoodRecoAnalyser() {
}

void LikelihoodRecoAnalyser::createHistograms() {
	histMan_->setCurrentHistogramFolder(histogramFolder_);

}

void LikelihoodRecoAnalyser::createTrees() {
	treeMan_->setCurrentFolder(histogramFolder_);
}

