/*
 * JetReader.h
 *
 *  Created on: Jun 25, 2010
 *      Author: lkreczko
 */

#ifndef JETREADER_H_
#define JETREADER_H_
#include "../Event.h"
#include "VariableReader.h"
#include "../RecoObjects/Jet.h"

namespace BAT {

class JetReader {
public:
	JetReader();
	JetReader(TChainPointer,
			JetAlgorithm::value algo = JetAlgorithm::Calo_AntiKT_Cone05);
	virtual ~JetReader();
	virtual const JetCollection& getJets(const bool isRealData);
	virtual void initialise();

protected:
	VariableReader<MultiDoublePointer> energyReader;
	VariableReader<MultiDoublePointer> JECReader;
	VariableReader<MultiDoublePointer> JECUncReader;
	VariableReader<MultiDoublePointer> L1OffJECReader;
	VariableReader<MultiDoublePointer> L2L3ResJECReader;
	VariableReader<MultiDoublePointer> L2RelJECReader;
	VariableReader<MultiDoublePointer> L3AbsJECReader;
	VariableReader<MultiDoublePointer> JERReader;
	VariableReader<MultiDoublePointer> JERSFReader;
	VariableReader<MultiDoublePointer> JERSFUpReader;
	VariableReader<MultiDoublePointer> JERSFDownReader;
	VariableReader<MultiDoublePointer> pxReader;
	VariableReader<MultiDoublePointer> pyReader;
	VariableReader<MultiDoublePointer> pzReader;
	VariableReader<MultiDoublePointer> energyRawReader;
	VariableReader<MultiDoublePointer> pxRawReader;
	VariableReader<MultiDoublePointer> pyRawReader;
	VariableReader<MultiDoublePointer> pzRawReader;
	VariableReader<MultiDoublePointer> massReader;
	VariableReader<MultiDoublePointer> chargeReader;

	VariableReader<MultiDoublePointer> matchedGeneratedJetEnergyReader;
	VariableReader<MultiDoublePointer> matchedGeneratedJetPxReader;
	VariableReader<MultiDoublePointer> matchedGeneratedJetPyReader;
	VariableReader<MultiDoublePointer> matchedGeneratedJetPzReader;

	VariableReader<MultiDoublePointer> matchedPartonEnergyReader;
	VariableReader<MultiDoublePointer> matchedPartonPxReader;
	VariableReader<MultiDoublePointer> matchedPartonPyReader;
	VariableReader<MultiDoublePointer> matchedPartonPzReader;
	VariableReader<MultiIntPointer> matchedPartonPdgIdReader;

	VariableReader<MultiDoublePointer> emfReader;
	VariableReader<unsigned int> numberOfJetsReader;
	VariableReader<MultiIntPointer> n90HitsReader;
	VariableReader<MultiDoublePointer> fHPDReader;

	VariableReader<MultiIntPointer> NODReader;
	VariableReader<MultiDoublePointer> CEFReader;
	VariableReader<MultiDoublePointer> NHFReader;
	VariableReader<MultiDoublePointer> NEFReader;
	VariableReader<MultiDoublePointer> CHFReader;
	VariableReader<MultiIntPointer> NCHReader;

	VariableReader<MultiDoublePointer> btagCSVv2Reader;
	VariableReader<MultiBoolPointer> isBtagReader;

	VariableReader<MultiDoublePointer> btagSFReader;
	VariableReader<MultiDoublePointer> btagSFUpReader;
	VariableReader<MultiDoublePointer> btagSFDownReader;

	VariableReader<MultiDoublePointer> tightbtagSFReader;
	VariableReader<MultiDoublePointer> tightbtagSFUpReader;
	VariableReader<MultiDoublePointer> tightbtagSFDownReader;

	VariableReader<MultiIntPointer> PartonFlavour;
	VariableReader<MultiIntPointer> HadronFlavour;

	VariableReader<MultiBoolPointer> isLooseReader;
	VariableReader<MultiBoolPointer> isMediumReader;
	VariableReader<MultiBoolPointer> isTightReader;

	JetCollection jets;
	JetAlgorithm::value usedAlgorithm;
	bool isRealData;
	virtual void readJets(bool);

};

}

#endif /* JETREADER_H_ */
