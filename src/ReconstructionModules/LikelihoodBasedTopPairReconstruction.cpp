/*
 * LikelihoodBasedTopPairReconstruction.cpp
 *
 *  Created on: 25 Aug 2011
 *      Author: kreczko
 */

#include "../../interface/ReconstructionModules/LikelihoodBasedTopPairReconstruction.h"
#include <boost/lexical_cast.hpp>

using namespace std;

namespace BAT {

LikelihoodBasedTopPairReconstruction::LikelihoodBasedTopPairReconstruction(const LeptonPointer lepton, const METPointer met, const JetCollection jets, const JetCollection bjets):
		solutions(),
		alreadyReconstructed(false),
		met(met),
		jets(jets),
		bjets(bjets),
		leptonFromW(lepton) {
}

bool LikelihoodBasedTopPairReconstruction::meetsInitialCriteria() const {
	return met != 0 && leptonFromW != 0 && jets.size() >= 4;
}

std::string LikelihoodBasedTopPairReconstruction::getDetailsOnFailure() const {
	std::string msg = "Initial Criteria not met: \n";
	if (leptonFromW == 0)
		msg += "Electron from W: not filled \n";
	else
		msg += "Electron from W: filled \n";

	if (met == 0)
		msg += "Missing transverse energy: not filled \n";
	else
		msg += "Missing transverse energy: filled \n";
	std::string nJets(boost::lexical_cast<std::string>(jets.size()));
	if (jets.size() < 4)
		msg += "Number of jets is too small:" + nJets + ", should be >= 4 \n";
	else
		msg += "Number of jets is OK:" + nJets + "\n";

	return msg;
}



LikelihoodBasedTopPairReconstruction::~LikelihoodBasedTopPairReconstruction() {
}

TtbarHypothesisCollection LikelihoodBasedTopPairReconstruction::getAllSolutions(){
	if(!alreadyReconstructed)
		reconstruct();

	SmallerDiscriminatorIsBetter solutionComparator = SmallerDiscriminatorIsBetter();
	std::sort(solutions.begin(), solutions.end(), solutionComparator);
	return solutions;
}

const TtbarHypothesisPointer LikelihoodBasedTopPairReconstruction::getBestSolution() {
	const TtbarHypothesisPointer bestSolution = getAllSolutions().front();//sorted by quality, front == best
	return bestSolution;
}

void LikelihoodBasedTopPairReconstruction::reconstruct() {
	typedef unsigned short ushort;

	// Loop b jets to get hadronic b candidate
	for (ushort hadBindex = 0; hadBindex < bjets.size(); ++hadBindex) {
		JetPointer hadBJet = bjets[hadBindex];
		if (!meetsHadronicBJetRequirement(hadBJet))
			continue;
		// Loop b jets to get leptonic b candidate
		for (ushort lepBindex = 0; lepBindex < bjets.size(); ++lepBindex) {
			JetPointer lepBJet = bjets[lepBindex];
			if (lepBindex == hadBindex || !meetsLeptonicBJetRequirement(lepBJet))
				continue;

			// Loop light jets to get jets from W
			for (ushort jet1Index = 0; jet1Index < jets.size(); ++jet1Index) {
				for (ushort jet2Index = 0; jet2Index < jets.size(); ++jet2Index) {
					JetPointer jet1FromW = jets[hadBindex];
					JetPointer jet2FromW = jets[hadBindex];

					if (jet2Index == jet1Index || !meetsJetFromWRequirement(jet1FromW, jet2FromW))
						continue;

					// Put all particles together
					TtbarHypothesisPointer solution(new TtbarHypothesis());
					//leptons
					solution->leptonFromW = leptonFromW;
					double neutrinoChi2 = -1;
					solution->neutrinoFromW = getNeutrinoSolution( neutrinoChi2 );
					//jets
					solution->hadronicBJet = hadBJet;
					solution->leptonicBjet = lepBJet;
					solution->jet1FromW = jet1FromW;
					solution->jet2FromW = jet2FromW;

					//combine reconstructed objects
					solution->combineReconstructedObjects();

					// Get discrimnant for this solution
					solution->discriminator = getDiscriminator(solution);

					// Store if event is physical (checks masses of tops and Ws > 0 )
					if (meetsGlobalRequirement(solution)){
						solutions.push_back(solution);
					}

					else
						continue;
				}
			}
		}
	}
	alreadyReconstructed = true;
}


ParticlePointer LikelihoodBasedTopPairReconstruction::getNeutrinoSolution( double& neutrinoChi2) {
	// Do the neutrino reconstruction
	return ParticlePointer( new Particle(1,1,1,1));
}

bool LikelihoodBasedTopPairReconstruction::meetsHadronicBJetRequirement(JetPointer hadBJet){
	return true;
}

bool LikelihoodBasedTopPairReconstruction::meetsLeptonicBJetRequirement(JetPointer lepBJet){
	return true;
}

bool LikelihoodBasedTopPairReconstruction::meetsJetFromWRequirement(JetPointer jet1, JetPointer jet2){
	return true;
}

bool LikelihoodBasedTopPairReconstruction::meetsGlobalRequirement(const TtbarHypothesisPointer solution){
	return solution->isPhysical();
}

double LikelihoodBasedTopPairReconstruction::getDiscriminator(const TtbarHypothesisPointer solution) const{
	// Given all the info for this solution/permutation, calculate the discriminant

	// Print csv of hadronic b jet and probablitity from correct b histogram
	JetPointer hadBJet = solution->hadronicBJet;
	double hadBJetCSV = hadBJet->getBTagDiscriminator(BAT::BtagAlgorithm::value::CombinedSecondaryVertexV2);

	double probCorrectB = Globals::csvCorrectPermHistogram->Interpolate( hadBJetCSV );
	cout << "CSV : " << hadBJetCSV << " prob : " << probCorrectB << endl;

	return 0;
}

} /* namespace BAT */
