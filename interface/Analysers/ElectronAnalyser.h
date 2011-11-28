/*
 * ElectronAnalyser.h
 *
 *  Created on: 3 Jul 2011
 *      Author: kreczko
 */

#ifndef ELECTRONANALYSER_H_
#define ELECTRONANALYSER_H_
#include "BasicAnalyser.h"

namespace BAT {

class ElectronAnalyser : public BasicAnalyser{
public:
    ElectronAnalyser(boost::shared_ptr<HistogramManager> histMan);
    virtual ~ElectronAnalyser();

    void analyse(const TopPairEventCandidate& ttbarEvent);
    void createHistograms();
};

}

#endif /* ELECTRONANALYSER_H_ */
