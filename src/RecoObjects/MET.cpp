/*
 * MET.cpp
 *
 *  Created on: Jun 25, 2010
 *      Author: lkreczko
 */

#include "../../interface/RecoObjects/MET.h"

namespace BAT {

MET::MET() :
		Particle(), //
		usedAlgorithm(METAlgorithm::MET), //
		significance_(-1), //
		sumET_(-1),
		metUncertaintyPt_() {
}

MET::MET(double ex, double ey) :
		Particle(sqrt(ex * ex + ey * ey), ex, ey, 0), //
		usedAlgorithm(METAlgorithm::MET), //
		significance_(-1), //
		sumET_(-1),
		metUncertaintyPt_() {

}

MET::~MET() {
}

bool MET::isGood() const {
	return et() > 20;
}

void MET::setUsedAlgorithm(METAlgorithm::value algo) {
	usedAlgorithm = algo;
}

METAlgorithm::value MET::getUsedAlgorithm() const {
	return usedAlgorithm;
}

void MET::setSignificance(double significance) {
	significance_ = significance;
}

void MET::setSumET(double sumET) {
	sumET_ = sumET;
}

double MET::significance() const {
	return significance_;
}

double MET::sumET() const {
	return sumET_;
}

double MET::resolution() const {
	return 0.5 * sqrt(sumET_);
}

bool MET::isMCOnlyMETType(unsigned int type) {
	// bool isMCOnlyMet = (type == METAlgorithm::patType1p2CorrectedPFMetJetResUp
	// 		|| type == METAlgorithm::patType1p2CorrectedPFMetJetResDown || type == METAlgorithm::GenMET);
	// isMCOnlyMet = isMCOnlyMet || type == METAlgorithm::patType1CorrectedPFMetJetResDown
	// 		|| type == METAlgorithm::patType1CorrectedPFMetJetResUp;
	// isMCOnlyMet = isMCOnlyMet || type == METAlgorithm::patPFMetJetResDown || type == METAlgorithm::patPFMetJetResUp;
	// return isMCOnlyMet;
	return false;
}

bool MET::isAvailableInNTupleVersion(unsigned int ntupleVersion, unsigned int type) {
// 	if (ntupleVersion > 8)
// 		return true;
// 	else if (ntupleVersion == 8)
// 		return type < METAlgorithm::recoMetPFlow;
// 	else if (ntupleVersion == 7)
// 		return type <= METAlgorithm::patType1p2CorrectedPFMetUnclusteredEnDown;
// 	else
// 		return type == METAlgorithm::patMETsPFlow || type == METAlgorithm::GenMET;
return true;
}

void MET::setMETUncertinaties( std::vector<double> newUncertainties ) {
	metUncertaintyPt_ = newUncertainties;
}

void MET::setMET_Px_Uncertinaties( std::vector<double> newUncertainties ) {
	metUncertaintyPx_ = newUncertainties;
}

void MET::setMET_Py_Uncertinaties( std::vector<double> newUncertainties ) {
	metUncertaintyPy_ = newUncertainties;
}

std::vector< double > MET::getAllMETUncertainties( ) {
	return metUncertaintyPt_;
}


METPointer MET::getMETForUncertainty( unsigned int unc_i ) {
	if ( unc_i < metUncertaintyPx_.size() ) {
		METPointer newMet = METPointer( new MET(metUncertaintyPx_.at( unc_i ), metUncertaintyPy_.at( unc_i )));
		return newMet;
	}
	else {
		return 0;
	}
}

}

